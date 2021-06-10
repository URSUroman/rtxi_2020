/*
	 Copyright (C) 2011 Georgia Institute of Technology, University of Utah, Weill Cornell Medical College

	 This program is free software: you can redistribute it and/or modify
	 it under the terms of the GNU General Public License as published by
	 the Free Software Foundation, either version 3 of the License, or
	 (at your option) any later version.

	 This program is distributed in the hope that it will be useful,
	 but WITHOUT ANY WARRANTY; without even the implied warranty of
	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	 GNU General Public License for more details.

	 You should have received a copy of the GNU General Public License
	 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <rtxi_config.h>
#include <cstring>
#include <string>
#include <unistd.h>
#include <compiler.h>
#include <debug.h>
#include <main_window.h>
#include <sstream>
#include <workspace.h>
#include <iostream>
#include <pthread.h>
#include <data_recorder.h>


//Added from .cpp file here to be able to add modifications from KZ
#include <math.h>
#include <fstream>
#include <sstream>



#define QFileExistsEvent            (QEvent::User+0)
#define QSetFileNameEditEvent       (QEvent::User+1)
#define QDisableGroupsEvent         (QEvent::User+2)
#define QEnableGroupsEvent          (QEvent::User+3)

#define TAG_SIZE 1024

struct param_hdf_t
{
    long long index;
    double value;
};

// Debug for event handling
QDebug operator<<(QDebug str, const QEvent * ev)
{
    static int eventEnumIndex = QEvent::staticMetaObject.indexOfEnumerator("Type");
    str << "QEvent";
    if (ev)
        {
            QString name = QEvent::staticMetaObject.enumerator(eventEnumIndex).valueToKey(ev->type());
            if (!name.isEmpty())
                str << name;
            else
                str << ev->type();
        }
    else
        str << (void*)ev;
    return str.maybeSpace();
}

struct find_daq_t
{
    int index;
    DAQ::Device *device;
};

static void findDAQDevice(DAQ::Device *dev,void *arg)
{
    struct find_daq_t *info = static_cast<struct find_daq_t *>(arg);
    if(!info->index)
        info->device = dev;
    info->index--;
}

namespace
{
void buildBlockPtrList(IO::Block *block, void *arg)
{
    std::vector<IO::Block *> *list = reinterpret_cast<std::vector<IO::Block *> *> (arg);
    list->push_back(block);
};

struct FileExistsEventData
{
    QString filename;
    int response;
    QWaitCondition done;
};

struct SetFileNameEditEventData
{
    QString filename;
    QWaitCondition done;
};

class InsertChannelEvent: public RT::Event
{
public:
    InsertChannelEvent(bool &, RT::List<DataRecorder::Channel> &,
                       RT::List<DataRecorder::Channel>::iterator, DataRecorder::Channel &);
    ~InsertChannelEvent(void);
    int callback(void);

private:
    bool &recording;
    RT::List<DataRecorder::Channel> &channels;
    RT::List<DataRecorder::Channel>::iterator end;
    DataRecorder::Channel &channel;
}; // class InsertChannelEvent

class RemoveChannelEvent: public RT::Event
{
public:
    RemoveChannelEvent(bool &, RT::List<DataRecorder::Channel> &,
                       DataRecorder::Channel &);
    ~RemoveChannelEvent(void);
    int callback(void);

private:
    bool &recording;
    RT::List<DataRecorder::Channel> &channels;
    DataRecorder::Channel &channel;
}; // class RemoveChannelEvent

class OpenFileEvent: public RT::Event
{
public:
    OpenFileEvent(QString &, AtomicFifo &);
    ~OpenFileEvent(void);
    int callback(void);

private:
    QString &filename;
    AtomicFifo &fifo;
}; // class OpenFileEvent

class StartRecordingEvent: public RT::Event
{
public:
    StartRecordingEvent(bool &, AtomicFifo &);
    ~StartRecordingEvent(void);
    int callback(void);

private:
    bool &recording;
    AtomicFifo &fifo;
}; // class StartRecordingEvent

class StopRecordingEvent: public RT::Event
{
public:
    StopRecordingEvent(bool &, AtomicFifo &);
    ~StopRecordingEvent(void);
    int callback(void);

private:
    bool &recording;
    AtomicFifo &fifo;
}; //class StopRecordingEvent

class AsyncDataEvent: public RT::Event
{
public:
    AsyncDataEvent(const double *, size_t, AtomicFifo &);
    ~AsyncDataEvent(void);
    int callback(void);

private:
    const double *data;
    size_t size;
    AtomicFifo &fifo;
}; // class AsyncDataEvent

class DoneEvent: public RT::Event
{
public:
    DoneEvent(AtomicFifo &);
    ~DoneEvent(void);
    int callback(void);

private:
    AtomicFifo &fifo;
}; // class DoneEvent
}; // namespace

InsertChannelEvent::InsertChannelEvent(bool &r, RT::List<DataRecorder::Channel> & l,
                                       RT::List<DataRecorder::Channel>::iterator e, DataRecorder::Channel &c) :
    recording(r), channels(l), end(e), channel(c)
{
}

InsertChannelEvent::~InsertChannelEvent(void)
{
}

int InsertChannelEvent::callback(void)
{
    if(recording)
        return -1;
    channels.insertRT(end, channel);
    return 0;
}

RemoveChannelEvent::RemoveChannelEvent(bool &r,	RT::List<DataRecorder::Channel> & l,
                                       DataRecorder::Channel &c) :
    recording(r), channels(l), channel(c)
{
}

RemoveChannelEvent::~RemoveChannelEvent(void)
{
}

int RemoveChannelEvent::callback(void)
{
    if (recording)
        return -1;
    channels.removeRT(channel);
    return 0;
}

OpenFileEvent::OpenFileEvent(QString &n, AtomicFifo &f) :
    filename(n), fifo(f)
{
}

OpenFileEvent::~OpenFileEvent(void)
{
}

int OpenFileEvent::callback(void)
{
    DataRecorder::data_token_t token;
    token.type = DataRecorder::OPEN;
    token.size = filename.length() + 1;
    token.time = RT::OS::getTime();

    //Added from KZ except second write: it was adapted as in original
    if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	    fifo.write(&token,sizeof(token));
	    fifo.write(filename.toLatin1().constData(), token.size);
    }

    //Commented from original
    // fifo.write(&token, sizeof(token));
    // fifo.write(filename.toLatin1().constData(), token.size);
    return 0;
}

StartRecordingEvent::StartRecordingEvent(bool &r, AtomicFifo &f) :
    recording(r), fifo(f)
{
}

StartRecordingEvent::~StartRecordingEvent(void)
{
}

int StartRecordingEvent::callback(void)
{
    DataRecorder::data_token_t token;
    recording = true;
    token.type = DataRecorder::START;
    token.size = 0;
    token.time = RT::OS::getTime();
    
    //Added from KZ
    if (!fifo.tooBig(sizeof(token)))
	    fifo.write(&token,sizeof(token));
    
    //Commented from original
    //fifo.write(&token, sizeof(token));
    return 0;
}

StopRecordingEvent::StopRecordingEvent(bool &r, AtomicFifo &f) :
    recording(r), fifo(f)
{
}

StopRecordingEvent::~StopRecordingEvent(void)
{
}

int StopRecordingEvent::callback(void)
{
    DataRecorder::data_token_t token;
    recording = false;
    token.type = DataRecorder::STOP;
    token.size = 0;
    token.time = RT::OS::getTime();

    //Added from KZ
    if (!fifo.tooBig(sizeof(token)))
      fifo.write(&token,sizeof(token));

    //Commented from original
    //fifo.write(&token, sizeof(token));
    return 0;
}

AsyncDataEvent::AsyncDataEvent(const double *d, size_t s, AtomicFifo &f) :
    data(d), size(s), fifo(f)
{
}

AsyncDataEvent::~AsyncDataEvent(void)
{
}

int AsyncDataEvent::callback(void)
{
    DataRecorder::data_token_t token;
    token.type = DataRecorder::ASYNC;
    token.size = size * sizeof(double);
    token.time = RT::OS::getTime();

    //Added from KZ
    if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
      fifo.write(&token,sizeof(token));
      fifo.write(data,token.size);
    }

    //Commented from original
    // fifo.write(&token, sizeof(token));
    // fifo.write(data, token.size);
    return 1;
}

DoneEvent::DoneEvent(AtomicFifo &f) :
    fifo(f)
{
}

DoneEvent::~DoneEvent(void)
{
}

int DoneEvent::callback(void)
{
    DataRecorder::data_token_t token;
    token.type = DataRecorder::DONE;
    token.size = 0;
    token.time = RT::OS::getTime();

    //Added from KZ
    if (!fifo.tooBig(sizeof(token)))
      fifo.write(&token,sizeof(token));

    //Commented from original
    //fifo.write(&token, sizeof(token));
    return 0;
}

DataRecorder::CustomEvent::CustomEvent(QEvent::Type type) : QEvent(type)
{
    data = 0;
}

void DataRecorder::CustomEvent::setData(void *ptr)
{
    data = ptr;
}

void * DataRecorder::CustomEvent::getData(void)
{
    return data;
}

void DataRecorder::startRecording(void)
{
    Event::Object event(Event::START_RECORDING_EVENT);
    if (RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

void DataRecorder::stopRecording(void)
{
    Event::Object event(Event::STOP_RECORDING_EVENT);
    if (RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

void DataRecorder::openFile(const QString& filename)
{
    Event::Object event(Event::OPEN_FILE_EVENT);
    event.setParam("filename", const_cast<char *> (filename.toLatin1().constData()));
    if (RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

void DataRecorder::postAsyncData(const double *data, size_t size)
{
    Event::Object event(Event::ASYNC_DATA_EVENT);
    event.setParam("data", const_cast<double *> (data));
    event.setParam("size", &size);
    if (RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

DataRecorder::Channel::Channel(void)
{
}

DataRecorder::Channel::~Channel(void)
{
}

//Added from KZ
static IO::channel_t chans[] = {
    {
        "Me",
        "",
        IO::OUTPUT,
    },
};
static size_t num_chans = sizeof(chans)/sizeof(chans[0]);


//Added from KZ: IO::Block("Connect",chans,num_chans) else in execute variable output(0) is not known
DataRecorder::Panel::Panel(QWidget *parent, size_t buffersize) :
  QWidget(parent), RT::Thread(RT::Thread::MinimumPriority), fifo(buffersize),  recording(false), IO::Block("Connect",chans,num_chans)
{
    setWhatsThis(
        "<p><b>Data Recorder:</b><br>The Data Recorder writes data to an HDF5 file format "
        "All available signals for saving to file are automatically detected. Currently "
        "loaded user modules are listed in the \"Block\" drop-down box. Available DAQ cards "
        "are listed here as /proc/analogy/devices. Use the \"Type\" and \"Channel\" drop-down boxes "
        "to select the signals that you want to save. Use the left and right arrow buttons to "
        "add these signals to the file. You may select a downsampling rate that is applied "
        "to the real-time period for execution (set in the System Control Panel). The real-time "
        "period and the data downsampling rate are both saved as metadata in the HDF5 file "
        "so that you can reconstruct your data correctly. The current recording status of "
        "the Data Recorder is shown at the bottom.</p>");


    // Make Mdi
    subWindow = new QMdiSubWindow;
    subWindow->setWindowIcon(QIcon("/usr/local/share/rtxi/RTXI-widget-icon.png"));
    subWindow->setAttribute(Qt::WA_DeleteOnClose);
    subWindow->setWindowFlags(Qt::CustomizeWindowHint | Qt::WindowCloseButtonHint |
                              Qt::WindowMinimizeButtonHint);
    MainWindow::getInstance()->createMdi(subWindow);

    // Create main layout
    QGridLayout *layout = new QGridLayout;

    // Create child widget and layout for channel selection
    channelGroup = new QGroupBox(tr("Channel Selection"));
    QVBoxLayout *channelLayout = new QVBoxLayout;

    // Create elements for channel box
    channelLayout->addWidget(new QLabel(tr("Block:")));
    blockList = new QComboBox;
    channelLayout->addWidget(blockList);
    QObject::connect(blockList,SIGNAL(activated(int)), this, SLOT(buildChannelList(void)));

    channelLayout->addWidget(new QLabel(tr("Type:")));
    typeList = new QComboBox;
    channelLayout->addWidget(typeList);
    typeList->addItem("Input");
    typeList->addItem("Output");
    typeList->addItem("Parameter");
    typeList->addItem("State");
    typeList->addItem("Event");
    QObject::connect(typeList,SIGNAL(activated(int)),this,SLOT(buildChannelList(void)));

    channelLayout->addWidget(new QLabel(tr("Channel:")));
    channelList = new QComboBox;
    channelLayout->addWidget(channelList);

    // Attach layout to child widget
    channelGroup->setLayout(channelLayout);

    // Create elements for arrow
    rButton = new QPushButton("Add");
    channelLayout->addWidget(rButton);
    QObject::connect(rButton,SIGNAL(released(void)),this,SLOT(insertChannel(void)));
    rButton->setEnabled(false);
    lButton = new QPushButton("Remove");
    channelLayout->addWidget(lButton);
    QObject::connect(lButton,SIGNAL(released(void)),this,SLOT(removeChannel(void)));
    lButton->setEnabled(false);

    // Timestamp
    stampGroup = new QGroupBox(tr("Tag Data"));
    QHBoxLayout *stampLayout = new QHBoxLayout;

    // Add timestamp elements
    timeStampEdit = new QLineEdit;
    stampLayout->addWidget(timeStampEdit);
    addTag = new QPushButton(tr("Tag"));
    stampLayout->addWidget(addTag);
    QObject::connect(addTag,SIGNAL(released(void)),this,SLOT(addNewTag(void)));

    // Attach layout to child widget
    stampGroup->setLayout(stampLayout);

    // Create child widget and layout
    sampleGroup = new QGroupBox(tr("Trial Metadata"));
    QHBoxLayout *sampleLayout = new QHBoxLayout;

    // create elements for sample box
    trialNumLbl = new QLabel;
    trialNumLbl->setText("Trial #:");
    trialNumLbl->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    sampleLayout->addWidget(trialNumLbl);
    trialNum = new QLabel;
    trialNum->setText("0");
    trialNum->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    sampleLayout->addWidget(trialNum);

    trialLengthLbl = new QLabel;
    trialLengthLbl->setText("Trial Length (s):");
    sampleLayout->addWidget(trialLengthLbl);
    trialLength = new QLabel;
    trialLength->setText("No data recorded");
    trialLength->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    sampleLayout->addWidget(trialLength);

    fileSizeLbl = new QLabel;
    fileSizeLbl->setText("File Size (MB):");
    sampleLayout->addWidget(fileSizeLbl);
    fileSize = new QLabel;
    fileSize->setText("No data recorded");
    sampleLayout->addWidget(fileSize);
    fileSize->setAlignment(Qt::AlignLeft | Qt::AlignVCenter);

    // Attach layout to child widget
    sampleGroup->setLayout(sampleLayout);

    // Create child widget and layout for file control
    fileGroup = new QGroupBox(tr("File Control"));
    QHBoxLayout *fileLayout = new QHBoxLayout;

    // Create elements for file control
    fileLayout->addWidget(new QLabel(tr("File Name:")));
    fileNameEdit = new QLineEdit;
    fileNameEdit->setReadOnly(true);
    fileLayout->addWidget(fileNameEdit);
    QPushButton *fileChangeButton = new QPushButton("Choose File");
    fileLayout->addWidget(fileChangeButton);
    QObject::connect(fileChangeButton,SIGNAL(released(void)),this,SLOT(changeDataFile(void)));

    fileLayout->addWidget(new QLabel(tr("Downsample \nRate:")));
    downsampleSpin = new QSpinBox(this);
    downsampleSpin->setMinimum(1);
    downsampleSpin->setMaximum(500);
    fileLayout->addWidget(downsampleSpin);
    QObject::connect(downsampleSpin,SIGNAL(valueChanged(int)),this,SLOT(updateDownsampleRate(int)));

    // Attach layout to child
    fileGroup->setLayout(fileLayout);

    // Create child widget and layout
    listGroup = new QGroupBox(tr("Currently Recording"));
    QGridLayout *listLayout = new QGridLayout;

    // Create elements for box
    selectionBox = new QListWidget;
    listLayout->addWidget(selectionBox,1,1,4,5);

    // Attach layout to child
    listGroup->setLayout(listLayout);

    // Creat child widget and layout for buttons
    buttonGroup = new QGroupBox;
    QHBoxLayout *buttonLayout = new QHBoxLayout;

    // Create elements for box
    startRecordButton = new QPushButton("Start Recording");
    QObject::connect(startRecordButton,SIGNAL(released(void)),this,SLOT(startRecordClicked(void)));
    buttonLayout->addWidget(startRecordButton);
    startRecordButton->setEnabled(false);
    stopRecordButton = new QPushButton("Stop Recording");
    QObject::connect(stopRecordButton,SIGNAL(released(void)),this,SLOT(stopRecordClicked(void)));
    buttonLayout->addWidget(stopRecordButton);
    stopRecordButton->setEnabled(false);
    closeButton = new QPushButton("Close");
    QObject::connect(closeButton,SIGNAL(released(void)),subWindow,SLOT(close()));
    buttonLayout->addWidget(closeButton);
    recordStatus = new QLabel;
    buttonLayout->addWidget(recordStatus);
    recordStatus->setText("Not ready.");
    recordStatus->setFrameStyle(QFrame::Panel | QFrame::Sunken);
    recordStatus->setAlignment(Qt::AlignHCenter | Qt::AlignVCenter);

    // Attach layout to group
    buttonGroup->setLayout(buttonLayout);

    // Attach child widgets to parent
    layout->addWidget(channelGroup, 0, 0, 1, 2);
    layout->addWidget(listGroup, 0, 2, 1, 4);
    layout->addWidget(stampGroup, 2, 0, 2, 6);
    layout->addWidget(fileGroup, 4, 0, 1, 6);
    layout->addWidget(sampleGroup, 5, 0, 1, 6);
    layout->addWidget(buttonGroup, 6, 0, 1, 6);

    setLayout(layout);
    setWindowTitle(QString::number(getID()) + " Data Recorder");

    // Set layout to Mdi
    subWindow->setWidget(this);
    subWindow->setFixedSize(subWindow->minimumSizeHint());
    show();

    // Register custom QEvents
    QEvent::registerEventType(QFileExistsEvent);
    QEvent::registerEventType(QSetFileNameEditEvent);
    QEvent::registerEventType(QDisableGroupsEvent);
    QEvent::registerEventType(QEnableGroupsEvent);

    // Build initial block list
    IO::Connector::getInstance()->foreachBlock(buildBlockPtrList, &blockPtrList);
    for (std::vector<IO::Block *>::const_iterator i = blockPtrList.begin(), end = blockPtrList.end(); i != end; ++i)
        blockList->addItem(QString::fromStdString((*i)->getName()) + " " + QString::number((*i)->getID()));

    // Setup thread sleep
    sleep.tv_sec = 0;
    sleep.tv_nsec = RT::System::getInstance()->getPeriod();

    // Check if FIFO is truly atomic for hardware arcstartecture
    if(!fifo.isLockFree())
        ERROR_MSG("DataRecorder::Panel: WARNING: Atomic FIFO is not lock free\n");

    // Build initial channel list
    buildChannelList();

    // Launch Recording Thread
    pthread_create(&thread, 0, bounce, this);
    counter = 0;
    downsample_rate = 1;
    prev_input = 0.0;
    count = 0;
    setActive(true);
}

// Destructor for Panel
DataRecorder::Panel::~Panel(void)
{
    Plugin::getInstance()->removeDataRecorderPanel(this);
    setActive(false);
    DoneEvent RTevent(fifo);
    while (RT::System::getInstance()->postEvent(&RTevent));
    pthread_join(thread, 0);
    for (RT::List<Channel>::iterator i = channels.begin(), end = channels.end(); i!= end;)
        delete &*(i++);
}



//************************
//Test: put start recording here to avoid :: in ispunct.It compiles. Not better behaviour in rtxi gui, it still crashes after 2 start/stop recordings
//***********************


// ***************************************************************
// **  **  **  **  **  **  My Code  **  **  **  **  **  **  **  **
// ***************************************************************

#include <../../../fftReal/ffft/FFTReal.h>	// Several of the FFTReal files have been modified
//#include <../electrode_resistance_measurement/electrode_resistance_measurement.h>
using namespace std;


// Parameters: Frequencies
const int SAMP_FREQ = 30.303 * 1000; 	// in Hertz
const int CUTOFF_FREQ = 6500;		// in Hz. The highest frequency considered for computing correlation coefficents
const int FREQ_LOW_CUT = 1000;		// in Hz. The lowest frequency considered.

// Parameters: RMS
const int RMS_PERIOD = 1000;		// in us
const double RMS_THRESHOLD = 0.15;//0.90;	// in Volts, threshold to start recording
const int MS_IN_RMS = 10;		// in ms, time length checked to start recording
const int RMS_TERM_LEN = 200;		// in ms, time length checked to stop recording
const double RMS_TERM_THRES = 0.07;//.85;	// in Volts, threshold to stop recording (if under, stops)

// Parameters: FFT
const int FFT_ARRAY_SIZE = 256; 	// Must be power of 2
const int FFT_PERIOD = 1000;		// in us

// Parameters: Recording: Lengths
const int PREDATA_LENGTH = 2;		// in s; Amount of Time Before Detection to Record. Also is seconds in Circular Buffer. Make >= 1 or program will crash.
const int RECORD_H5_SECS = 5;		// in s; Min Num Seconds/ Trial Recorded
const int MAX_FILE_SECS = 25;		// in s; Max sec.s recorded, not including predata

// Parameters: Spectrogram/ Template
const char TEMP_FILE[45] = "Template_syl_a_B2.txt"; //"LowerHalfTemp.Det";	// Template File (Raw Data)
const int MS_IN_SPECT = 55;//55;		// in ms
const int SPECT_PERIOD = 1000;		// in us

// Parameters: Trigger
double corr_coeff_trig = 0.7;//0.74;		// Correlation to Template needed for Detection

// Parameters: Output
const int OUTPUT_LENGTH = 100; //32;		// in ms; output file length
const char PLAYBACK_FILE[45] = "Noise.txt"; //"thisSyl"; // File for Feedback
const int DEAD_TIME = 120;//50; 		// in ms. Default: OUTPUT_LENGTH + 20

// End Parameters                      ******                      ******
// End Parameters                      ******                      ******
// End Parameters                      ******                      ******

// Circular Buffer
const int CB_ARRAY_SIZE = PREDATA_LENGTH * SAMP_FREQ;
double circularBuffer[2][CB_ARRAY_SIZE];
double circularBuffer_sqrd[CB_ARRAY_SIZE];
double circularBuffer_flt[CB_ARRAY_SIZE];
double circularBuffer_bp_flt[CB_ARRAY_SIZE];
double predata[2][CB_ARRAY_SIZE];
double predata_flt[CB_ARRAY_SIZE];

// RMS Range
const int POINTS_IN_RMS = SAMP_FREQ / 1000 * MS_IN_RMS;
const int POINTS_IN_RMS_TERM = SAMP_FREQ / 1000 * RMS_TERM_LEN;

// FFT Objects
const long FFT_LEN = long(FFT_ARRAY_SIZE);	// Documentation in fftReal Folder
ffft::FFTReal <double> fft_object (FFT_LEN);
const int CUTOFF_INDEX = floor(CUTOFF_FREQ * FFT_ARRAY_SIZE / SAMP_FREQ);
const int BEGIN_INDEX = ceil(FREQ_LOW_CUT * FFT_ARRAY_SIZE / SAMP_FREQ);
double window[FFT_ARRAY_SIZE];

// Spectrogram Objects
const int SPECT_ARRAY_SIZE = MS_IN_SPECT * 1000 / FFT_PERIOD;
double spectBuffer[SPECT_ARRAY_SIZE][FFT_ARRAY_SIZE/2];
double spectMean[SPECT_ARRAY_SIZE];
double temp_spect[SPECT_ARRAY_SIZE][FFT_ARRAY_SIZE/2];
double tempLen;

// Output Objects
const int OUT_LEN_PTS = SAMP_FREQ / 1000 * OUTPUT_LENGTH;
double playback[OUT_LEN_PTS];

// Record Objects
const int RECORD_MAX = SAMP_FREQ * RECORD_H5_SECS;
const int MFS_PTS = SAMP_FREQ * MAX_FILE_SECS;
double dataStorage[2][MFS_PTS];
double dataStorage_flt[MFS_PTS];
int trialNum = 0;

// Counters
int cBcounter = 0;
int record_counter = 0;
int sBcounter = 0;
int pBcounter = 0;
int dScounter = 0;
int deadCounter = 0;

// Booleans
bool initialized = false;
bool isRunning = false;
bool rmsHit = false;
bool zapEm = false;
bool writeData = false;
bool isInteresting = true;
bool Template_found = false;
bool pre_zapEm = false;

// Physical Constants
const double PI = 3.1415926535897932;
char cwd[PATH_MAX];

//Additional variables/constants with respect to KZ
//Time variable for labeling file names
int write_control;
ofstream file_to_write;
ofstream file_to_write_flt;
time_t current_time;
string file_path = "../../Documents/Recordings/Song/Raw_Recordings/";//"../../Documents/Recordings/Song/";
string file_path_flt = "../../Documents/Recordings/Song/Raw_Recordings_flt/";
string file_name;
string file_name_flt;
//const char* file_name;
stringstream ss;
stringstream ss_flt;
string file_id;
string file_id_flt;

int max_h5_counter = SAMP_FREQ;
int h5_counter = 0;
double max_corr_coef = 0;

//***********************************************
//Delay for hte noise output
//***********************************************
int Dt_ms = 1;
int ms_counter = 0;
int max_ms_counter = Dt_ms*30; //Roughly the number of samples in Dt_ms ms
int min_dur_sil_ms = 1;
int min_dur_sil = 30*min_dur_sil_ms;
int min_dur_sil_count = 0;
//End addintional variables

//Filter parameters
const char FILTER_FILE_H[45] = "Filter_h.txt"; //Values for filling the filter
const char FILTER_FILE_B[45] = "Filter_b.txt"; //Values for filling the filter
const int FILTER_SIZE_H = 303;
const int FILTER_SIZE_B = 513;
double filter_h[FILTER_SIZE_H]={0};
double filter_b[FILTER_SIZE_B]={0};
double filter_sz_b = double(FILTER_SIZE_B);
double filter_sz_h = double(FILTER_SIZE_H);
double data_flt = 0;
double data_sqrd = 0;
//********************************************
//Amplitude threshold for syllable duration
//********************************************
double Amp_th = 1.6e-4;
double time_counter = 0;
//*******************************************
//Threshold for syllable duration
//*******************************************
double syl_th_duration = 1000000;//corresponds to the nb of samples at 30303 Hz.
const int Nb_syllables = 200;
double syllables_durations[Nb_syllables] = {0};
int syl_counter = 0;


// Functions
void Initialize();
void FftTemp(double temp_data[], int, int&, double&);
void Command_Code(double, double);
void Begin();
double RmsValue(int);
void FftGo();
void Compare_Spects();
void StopStart();
void AppendMe(double data[]);
void Reset();

// ***********************************************************
// ********************* Start Functions *********************
// ***********************************************************

// Initialize All The Things!
void Initialize()
{
	// Initialize Window
	for (int i = 0; i != FFT_ARRAY_SIZE; i++)
		window[i] = (-cos(2*PI*i/FFT_ARRAY_SIZE) + 1);

	// Declare Array and Counters
	const int TD_SIZE = FFT_ARRAY_SIZE + int(SAMP_FREQ / 1000) * MS_IN_SPECT;
	double temp_data[TD_SIZE];
	int tDcounter = 0;
	int tScounter = 0;

	//No write
	write_control = 0;

	// // Get File Data
	// ifstream playbackFile;
	// playbackFile.open(PLAYBACK_FILE);
	// for (int i = 0; i != OUT_LEN_PTS; i++)
	// 	playbackFile >> playback[i];
	// playbackFile.close();

	ifstream playbackFile;
	playbackFile.open(PLAYBACK_FILE);
	if(playbackFile.is_open()){
	  for (int i = 0; i != OUT_LEN_PTS; i++)
	    playbackFile >> playback[i];
	}
	else{
	  cout<<"Unable to open playback file. File open is: "<<playbackFile.is_open()<<endl;
	  if (getcwd(cwd, sizeof(cwd)) != NULL) {
	    cout<<"Current working pback dir: "<<cwd<<endl;
	  }
	}
	playbackFile.close();


	// // Get Temp Data
	// ifstream tempFile;
	// tempFile.open(TEMP_FILE);
	// for (int i = 0; i != TD_SIZE; i++) 
	// 	tempFile >> temp_data[i];
	// tempFile.close();


	ifstream tempFile;
	tempFile.open(TEMP_FILE);
	if(tempFile.is_open()){
	  for (int i = 0; i != TD_SIZE; i++) 
	    tempFile >> temp_data[i];
	}
	else{
	  cout<<"Unable to open template file. File open is: "<<tempFile.is_open()<<endl;
	  if (getcwd(cwd, sizeof(cwd)) != NULL) {
	    cout<<"Current working tpl dir: "<<cwd<<endl;
	  }
	}
	tempFile.close();


	//Fill filter
	ifstream filterFile;
	filterFile.open(FILTER_FILE_B);
	if(filterFile.is_open()){
	  for (int i = 0; i != FILTER_SIZE_B; i++)
	    filterFile >> filter_b[i];
	}
	else{
	  cout<<"Unable to open filter file. File open is: "<<filterFile.is_open()<<endl;
	  if (getcwd(cwd, sizeof(cwd)) != NULL) {
	    cout<<"Current working filter dir: "<<cwd<<endl;
	  }
	}
	filterFile.close();

	//Fill filter
	filterFile.open(FILTER_FILE_H);
	if(filterFile.is_open()){
	  for (int i = 0; i != FILTER_SIZE_H; i++)
	    filterFile >> filter_h[i];
	}
	else{
	  cout<<"Unable to open filter file. File open is: "<<filterFile.is_open()<<endl;
	  if (getcwd(cwd, sizeof(cwd)) != NULL) {
	    cout<<"Current working filter dir: "<<cwd<<endl;
	  }
	}
	filterFile.close();



	// Fill Spectrogram, Do FFT's and Calculate Mean
	double columnMean;
	double mean = 0;
	for (tDcounter = FFT_ARRAY_SIZE; tDcounter != TD_SIZE; tDcounter++)
		if (( (tDcounter - FFT_ARRAY_SIZE) % (SAMP_FREQ / (1000000 / FFT_PERIOD) ) == 0)) {
			FftTemp(temp_data, tDcounter, tScounter, columnMean);
			mean = mean + columnMean;
		}

	// Subtract the Mean
	mean = mean / SPECT_ARRAY_SIZE;
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		for (int j = BEGIN_INDEX;  (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			temp_spect[i][j] = temp_spect[i][j] - mean;

	// Calculate Template Length
	tempLen = 0;
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		for (int j = BEGIN_INDEX;   (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			tempLen = tempLen + pow(temp_spect[i][j],2);
	tempLen = sqrt(tempLen);

	// Normalize Template
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		for (int j = BEGIN_INDEX;   (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			temp_spect[i][j] = temp_spect[i][j] / tempLen;

	initialized = true;

	cout << "\nSpectral Analysis Engine Ready.\n\n"; // open RTXI from the terminal to see cout messages

}

// FFT the Template
void FftTemp(double temp_data[], int tDcounter, int& tScounter, double& mean)
{
	// Fill FFT-In Array
	double fftBuffer[FFT_ARRAY_SIZE];
	for (int i = 0; i != FFT_ARRAY_SIZE; i++)
		fftBuffer[i] = temp_data[tDcounter - (FFT_ARRAY_SIZE - i)];

	// Window FFT-In
	for (int i = 0; i != FFT_ARRAY_SIZE; i++)
		fftBuffer[i] = fftBuffer[i] * window[i];

	// Do FFT [Documentation included in FFTReal]
	double transformedFFT[FFT_ARRAY_SIZE];
	fft_object.do_fft (transformedFFT, fftBuffer);

	// Calculate Log|FFT'd|
	const int HALF_FFT = FFT_ARRAY_SIZE / 2;
	transformedFFT[0] = abs(transformedFFT[0]);
	for (int i = BEGIN_INDEX; (i != HALF_FFT) && (i != CUTOFF_INDEX); i++)
		transformedFFT[i] = log( sqrt(pow(transformedFFT[i],2) + pow(transformedFFT[i+HALF_FFT],2) ) );

	// Calculate Mean of FFT'd
	mean = 0;
	for (int i = BEGIN_INDEX; (i != HALF_FFT) && (i != CUTOFF_INDEX); i++)
		mean = mean + transformedFFT[i];
	mean = mean / HALF_FFT;

	// Fill Spectrogram Array
	for (int i = BEGIN_INDEX; (i != HALF_FFT) && (i != CUTOFF_INDEX); i++) 
		temp_spect[tScounter][i] = transformedFFT[i];
	tScounter++;
}

// Main Code Block
void Command_Code(double dataCh1, double dataCh2)
{
	// Begin Analysis
	if (!isRunning) Begin();

	// Fill Buffer
	circularBuffer[0][cBcounter] = dataCh1;
	circularBuffer[1][cBcounter] = dataCh2;


        //Band pass filter with freq cutoffs 1000;8000
        data_flt = 0;
	if ((cBcounter >= (FILTER_SIZE_B-1))) {
	  for (int i = 0; i != FILTER_SIZE_B; i++)
	    data_flt = data_flt + circularBuffer[0][cBcounter-i]*filter_b[i];
	    //filtered_point = filtered_point + circularBuffer_sqrd[cBcounter-i]*filter[i];
	}
	else {
	  for (int i = 0; i != (cBcounter+1); i++)
	    data_flt = data_flt + circularBuffer[0][cBcounter-i]*filter_b[i];
	    //filtered_point = filtered_point + circularBuffer_sqrd[cBcounter-i]*filter[i];
	  for (int i = 0; i != (FILTER_SIZE_B-cBcounter); i++)
	    data_flt = data_flt + circularBuffer[0][CB_ARRAY_SIZE-1-i]*filter_b[cBcounter+1+i];
	    //filtered_point = filtered_point + circularBuffer_sqrd[CB_ARRAY_SIZE-1-i]*filter[cBcounter+i];
	}
        //circularBuffer_filtd[cBcounter] = filtered_point;
	circularBuffer_bp_flt[cBcounter] = data_flt;


	//Square the bp filtered signal
	circularBuffer_sqrd[cBcounter] = circularBuffer_bp_flt[cBcounter]*circularBuffer_bp_flt[cBcounter];

	//Filter the squared signal
	// if ((cBcounter >= FILTER_SIZE)) {
	//   data_flt = data_flt/filter_sz + (circularBuffer_sqrd[cBcounter]-circularBuffer_sqrd[cBcounter-FILTER_SIZE]);
	// }
	// else {
	//   data_flt = data_flt/filter_sz + (circularBuffer_sqrd[cBcounter]-circularBuffer_sqrd[CB_ARRAY_SIZE - FILTER_SIZE + cBcounter]);
	// }

	//Filter the squared signal
	// if ((cBcounter >= FILTER_SIZE)) {
	//   data_flt = data_flt + (circularBuffer_sqrd[cBcounter] - circularBuffer_sqrd[cBcounter-FILTER_SIZE])/filter_sz;
	// }
	// else {
	//   data_flt = data_flt + (circularBuffer_sqrd[cBcounter] - circularBuffer_sqrd[CB_ARRAY_SIZE - FILTER_SIZE + cBcounter])/filter_sz;
	// }
	// circularBuffer_flt[cBcounter] = data_flt;



        //High pass filter
        data_flt = 0;
	if ((cBcounter >= (FILTER_SIZE_H-1))) {
	  for (int i = 0; i != FILTER_SIZE_H; i++)
	    data_flt = data_flt + circularBuffer_sqrd[cBcounter-i]*filter_h[i];
	    //filtered_point = filtered_point + circularBuffer_sqrd[cBcounter-i]*filter[i];
	}
	else {
	  for (int i = 0; i != (cBcounter+1); i++)
	    data_flt = data_flt + circularBuffer_sqrd[cBcounter-i]*filter_h[i];
	    //filtered_point = filtered_point + circularBuffer_sqrd[cBcounter-i]*filter[i];
	  for (int i = 0; i != (FILTER_SIZE_H-cBcounter); i++)
	    data_flt = data_flt + circularBuffer_sqrd[CB_ARRAY_SIZE-1-i]*filter_h[cBcounter+1+i];
	    //filtered_point = filtered_point + circularBuffer_sqrd[CB_ARRAY_SIZE-1-i]*filter[cBcounter+i];
	}
        //circularBuffer_filtd[cBcounter] = filtered_point;
	circularBuffer_flt[cBcounter] = data_flt;
	



	//Check whether the amplitude is above Amp_th.
	if(( circularBuffer_flt[cBcounter]> Amp_th))
	  {
	    if(min_dur_sil_count>0){//amplitude exceded thrsh before silence lasted a minimum duration. Add min_dur_sil_count to time_counter
	      time_counter = time_counter + min_dur_sil_count;
	      min_dur_sil_count = 0;
	    }
	    else{//amplitude exceded thrsh after silence lasted a minimum duration
	      time_counter++;
	    }
	  }
	else
	  {
	    if(time_counter > 0) //A syllable has just finished. time_counter contains the duration of the syllable
	      {
		//Check that the silence has a minimum duration of min_dur_sil
		if(min_dur_sil_count<min_dur_sil)
		  {
		    min_dur_sil_count++;
		  }
		else{//The silence lasted longer that thrshd duration
		  if(Template_found) //The syllable was the right one (  && (zapEm == false) && (pre_zapEm == false))
		    {
		      //Record the syllable duration
		      //syllables_durations[syl_counter] = time_counter;
		      //syl_counter++;
		      cout<<"tim_cntr(syl_dur): "<<time_counter<<endl;
		      if((time_counter < syl_th_duration)) //The syllable is shorther than a threshold
			{
			  //pre_zapEm = true; //Allow noise output with delay
			  zapEm = true;
			}
		      Template_found = false;
		    }
		  time_counter = 0;
		  min_dur_sil_count=0;
		}
	      }
	  }

	
	// Do RMS (when not recording) and Grab Pre-Data
	if (!rmsHit)
		if ((  cBcounter % (SAMP_FREQ / (1000000 / RMS_PERIOD)) == ( (SAMP_FREQ / (1000000 / RMS_PERIOD)) - 1) ))
			if ((RmsValue(POINTS_IN_RMS) >= RMS_THRESHOLD)) {
				rmsHit = true;
				for (int i = 0, j = cBcounter+1; i <  CB_ARRAY_SIZE; i++, j++) {
					predata[0][i] = circularBuffer[0][j];
					predata[1][i] = circularBuffer[1][j];
					//Filtered data
					predata_flt[i] = circularBuffer_flt[j];
					if (j == CB_ARRAY_SIZE - 1) j = -1;
				}
			}

	// Do FFT
	if (rmsHit){
	  if (((cBcounter+1) % (SAMP_FREQ / (1000000 / FFT_PERIOD)) == ( (SAMP_FREQ / (1000000 / FFT_PERIOD)) - 1))){
	    FftGo();
	  }
	}

	// Compare Spectrogram with Template
	if (rmsHit){
	  if (((cBcounter+2) % (SAMP_FREQ / (1000000 / SPECT_PERIOD)) == ( (SAMP_FREQ / (1000000 / SPECT_PERIOD)) - 1))){
	    Compare_Spects();
	  }
	}
	// Record H5 for X seconds
	if (rmsHit) {
		if ((record_counter++ == RECORD_MAX)) {
			if ((RmsValue(POINTS_IN_RMS_TERM) < RMS_TERM_THRES)) {
				rmsHit = false;
				record_counter = 0;
				writeData = true;
			} else record_counter = RECORD_MAX - SAMP_FREQ;
		}
	}

	// Make Buffer wrap around
	if ((cBcounter == CB_ARRAY_SIZE - 1)) {
		cBcounter = 0;
	} else cBcounter++;
}

// Begin Analysis
void Begin()
{
	isRunning = true;

	cout << "\nBegin Analysis.\n\n"; // open RTXI from the terminal to see cout messages
}

// Compute RMS Value
double RmsValue(int length) 
{
	double RMS = 0;
	int rmsLoop;

	// No wrap around
	if ((cBcounter >= length))
		for (rmsLoop = cBcounter- length; rmsLoop != cBcounter; rmsLoop++)
			RMS = RMS + pow(circularBuffer[0][rmsLoop],2);

	// Wrap Around Buffer
	else {
		int rmsFiller = 0;
		for (rmsLoop = CB_ARRAY_SIZE - (length - cBcounter); rmsLoop != CB_ARRAY_SIZE; rmsLoop++, rmsFiller++)
			RMS = RMS + pow(circularBuffer[0][rmsLoop],2);
		for (rmsLoop = 0; rmsFiller != length; rmsFiller++, rmsLoop++)
			RMS = RMS + pow(circularBuffer[0][rmsLoop],2);
	}

	// Finalize Value
	//RMS = sqrt(RMS/ (SAMP_FREQ / (1000000 / RMS_PERIOD)) );
	RMS = sqrt(RMS/length);
	//cout<<"RMS: "<<RMS<<"\n";
	return RMS;
}

// Do FFT
void FftGo()
{
	// Fill FFT-In Array 
	int k = 0;
	double fftBuffer[FFT_ARRAY_SIZE];
	if ((cBcounter >= FFT_ARRAY_SIZE)) {
		for (int i = cBcounter - FFT_ARRAY_SIZE; i != cBcounter; i++, k++)
			fftBuffer[k] = circularBuffer[0][i];
	}
	else {
		for (int i = CB_ARRAY_SIZE - (FFT_ARRAY_SIZE - cBcounter); i != CB_ARRAY_SIZE; i++, k++)
			fftBuffer[k] = circularBuffer[0][i];
		for (int i = 0; k != FFT_ARRAY_SIZE; i++, k++)
			fftBuffer[k] = circularBuffer[0][i];
	}

	// Window FFT-In 
	for (int i = 0; i != FFT_ARRAY_SIZE; i++)
		fftBuffer[i] = fftBuffer[i] * window[i];

	// Do FFT [Full documentation in FFTReal Folder]
	double transformedFFT[FFT_ARRAY_SIZE];
	fft_object.do_fft (transformedFFT, fftBuffer);

	// Calculate Log|FFT'd| 
	const int HALF_FFT = FFT_ARRAY_SIZE / 2;
	transformedFFT[0] = abs(transformedFFT[0]);
	for (int i = BEGIN_INDEX; (i != HALF_FFT) && (i != CUTOFF_INDEX); i++)
		transformedFFT[i] = log( sqrt(pow(transformedFFT[i],2) + pow(transformedFFT[i+HALF_FFT],2) ) );

	// Calculate Mean of FFT'd 
	spectMean[sBcounter] = 0;
	for (int i = BEGIN_INDEX; (i != HALF_FFT) && (i != CUTOFF_INDEX); i++)
		spectMean[sBcounter] = spectMean[sBcounter] + transformedFFT[i];
	spectMean[sBcounter] = spectMean[sBcounter] / HALF_FFT;

	// Fill Spectrogram Array 
	for (int i = BEGIN_INDEX; (i != HALF_FFT) && (i != CUTOFF_INDEX); i++)
		spectBuffer[sBcounter][i] = transformedFFT[i];
	if ((sBcounter != SPECT_ARRAY_SIZE - 1)) sBcounter++;
	else sBcounter = 0;
}

// Compare Current Spectrogram with Template
void Compare_Spects()
{
	// Create Current Spectrogram 
	double currentSpect[SPECT_ARRAY_SIZE][FFT_ARRAY_SIZE / 2];
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		for (int j = BEGIN_INDEX; (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			currentSpect[i][j] = spectBuffer[i][j];

	// Subtract Mean from Current 
	double mean = 0;
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		mean = mean + spectMean[i];
	mean = mean / SPECT_ARRAY_SIZE;
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		for (int j = BEGIN_INDEX; (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			currentSpect[i][j] = currentSpect[i][j] - mean;

	// Calculate Length of Current 
	double spectLen = 0;
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++)
		for (int j = BEGIN_INDEX; (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			spectLen = spectLen + pow(currentSpect[i][j],2);
	spectLen = sqrt(spectLen);

	// Dot Normalized Current with Template 
	int k = 0;
	double correlation_coefficient = 0;
	for (int i = sBcounter; i != SPECT_ARRAY_SIZE; i++, k++)
		for (int j = BEGIN_INDEX; (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			correlation_coefficient = correlation_coefficient + currentSpect[i][j] * temp_spect[k][j];
	for (int i = 0; k != SPECT_ARRAY_SIZE; i++, k++)
		for (int j = BEGIN_INDEX; (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			correlation_coefficient = correlation_coefficient + currentSpect[i][j] * temp_spect[k][j];

	// Normalize 
	correlation_coefficient = correlation_coefficient / spectLen;
	//max_corr_coef = max(correlation_coefficient,max_corr_coef);
	//cout<<"correlation_coefficient: "<<correlation_coefficient<<endl;
	//cout<<"currentSpect: "<<spectLen<<endl;
	

	// Zap Em? (output)
	if (correlation_coefficient >= corr_coeff_trig && deadCounter == 0) {
		isInteresting = true;
		Template_found = true;
		//zapEm = true;
	}
}

// On to the Next Trial!
void StopStart()
{
	DataRecorder::stopRecording();
	DataRecorder::startRecording();
	isInteresting = false;
}

// Append Me
void DataRecorder::Panel::AppendMe(double thisData[])
{
	H5PTappend(file.cdata,1,thisData);
}

// Reset Everything
void Reset()
{
	// Counters
	cBcounter = 0;
	record_counter = 0;
	sBcounter = 0;
	pBcounter = 0;
	dScounter = 0;
	deadCounter = 0;

	// Booleans
	isRunning = false;
	rmsHit = false;
	zapEm = false;
	writeData = false;
	isInteresting = true;
	Template_found = false;
	pre_zapEm = false;

	// Arrays
	for (int i = 0; i != SPECT_ARRAY_SIZE; i++) {
		spectMean[i] = 0;
		for (int j = BEGIN_INDEX;  (j != FFT_ARRAY_SIZE / 2) && (j != CUTOFF_INDEX); j++)
			spectBuffer[i][j] = 0;
	}
	for (int i = 0; i != CB_ARRAY_SIZE; i++) {
		circularBuffer[0][i] = 0;
	}

	cout << "Analysis Complete.\n\n"; // open RTXI from the terminal to see cout messages
}

// ***************************************************************
// **  **  **  **  **  End My_Code  **  **  **  **  **  **  **  **
// ***************************************************************


// Note: Execute has also been heavily modified, and there are other spotted changes throughout this file, fifo.cpp, fifo.h and data_recorder.h

void DataRecorder::Panel::execute(void) { 
   if(recording && !counter++) {

	// Get data Pts
        double data[channels.size()];
        size_t n = 0;

	data_token_t token;
	//token.type = SYNC;
	token.size = channels.size() * sizeof(double);
	
        for(RT::List<Channel>::iterator i = channels.begin(),end = channels.end();i != end;++i)
            if(i->block) data[n++] = i->block->getValue(i->type,i->index);

	// Activate Command Code
	if (!writeData) Command_Code(data[0],data[1]);

	//In case we need to delay noise signal, use pre_zapEm
	// if(pre_zapEm && (ms_counter < max_ms_counter)){
	//   //cout<<"max_corr_coeff: "<<max_corr_coef<<endl;
	//   ms_counter++;
	// }
	// if(ms_counter==max_ms_counter){
	//   //zapEm = true;
	//   pre_zapEm = false;
	// }


	
	// Output Playback File
	if (zapEm) {
		if (pBcounter < OUT_LEN_PTS)
		  {
		    output(0) = playback[pBcounter];
		    //cout<<"output playback"<<endl;
		  }
		else if (pBcounter == OUT_LEN_PTS)
			output(0) = 0.0;
		if ((++pBcounter == OUT_LEN_PTS+1)) {
			pBcounter = 0;
			zapEm = false;
			deadCounter = 1;
			ms_counter=0;
		}
	}

	// Write Data
	if (writeData) {
		// if (isInteresting) {
		// 	double appendData[2];
		// 	for (int i = 0; i != CB_ARRAY_SIZE; i++) {
		// 		appendData[0] = predata[0][i];
		// 		appendData[1] = predata[1][i];
		// 		DataRecorder::Panel::AppendMe(appendData);
		// 	}
		// 	for (int i = 0; i != dScounter; i++) {
		// 		appendData[0] = dataStorage[0][i];
		// 		appendData[1] = dataStorage[1][i];
		// 		DataRecorder::Panel::AppendMe(appendData);
		// 	}
		// 	cout << "  Trial " << trialNum << " Complete.\n\n"; // open RTXI from the terminal to see cout messages
		// 	//StopStart();
		// }

	  //Write data to fifos to be written to files in the low priority thread
	  //***********************************************************************************
	  //                  Coment section bellow to not write data to files
	  //***********************************************************************************
	  if (isInteresting) {
	    //double appendData[2 * channels.size()];
	    double appendData[channels.size()];

	    //Put _START label in fifo
	    token.type = _START;
	    //token.size = 2 * channels.size() * sizeof(double); //Data point and Data_flt point are doubles . channels.size is normally 1, that's why we have 2*channels.size
	    token.size = channels.size() * sizeof(double);
	    fifo.write(&token, sizeof(token));

	    //Put data into fifo
	    token.type = SYNC;
	    for (int i = 0; i != (CB_ARRAY_SIZE-1); i++) {
	      //appendData[0] = predata[0][i];
	      //appendData[1] = predata[1][i];
	      appendData[0] = predata[0][i];
	      //appendData[1] = predata_flt[i];
	      if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	  	fifo.write(&token, sizeof(token));
	  	fifo.write(appendData, sizeof(appendData));
	      }
	      // else
	      //   {
	      //     cout<<"fifo.tooBig for predata: "<<endl;
	      //   }
	    }
			
	    for (int i = 0; i != dScounter; i++) {
	      //appendData[0] = dataStorage[0][i];
	      //appendData[1] = dataStorage[1][i];
	      appendData[0] = dataStorage[0][i];
	      //appendData[1] = dataStorage_flt[i];
	      if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	  	fifo.write(&token, sizeof(token));
	  	fifo.write(appendData, sizeof(appendData));
	      }
	      // else
	      //   {
	      //     cout<<"fifo.tooBig for dataStorage: "<<endl;
	      //   }

	    }
	   cout << "  Trial " << trialNum << " Complete.\n\n"; // open RTXI from the terminal to see cout messages

	    //Put _STOP token in fifo
	    token.type = _STOP;
	    fifo.write(&token, sizeof(token));

	    write_control = 1;
			
	    //StopStart();
	  }
	  //***********************************************************************************
	  //                  Coment section above to not write data to files
	  //***********************************************************************************

	  writeData = false;
	  dScounter = 0;
	}

	
	// Store Data to be Written (several seconds) Later
	if (rmsHit && (dScounter < MFS_PTS)) {
		dataStorage[0][dScounter] = data[0];
		dataStorage[1][dScounter] = data[1]; //dScounter++ was here
		//Filtered data
		dataStorage_flt[dScounter++] = circularBuffer_flt[cBcounter-1];
	}
    }
    // OP Code
    count++;
    counter %= downsample_rate;

    // Reset, Ground Output
    if ((!recording && isRunning)) {
	Reset();
	output(0) = 0.0;
    }
    
    // Intialize Everything
    if (!initialized){
      Initialize();
      cout<<"Fifo size: "<<fifo.GetFifoSize()<<endl;
    }

    // Dead Time
    if (deadCounter != 0) deadCounter++;
    if (deadCounter >= DEAD_TIME * (SAMP_FREQ/ 1000)) deadCounter = 0;
}



//Commented from original
// Execute loop (original)
/*
void DataRecorder::Panel::execute(void)
{
    if (recording && !counter++)
        {
            data_token_t token;
            double data[channels.size()];

            size_t n = 0;
            token.type = SYNC;
            token.size = channels.size() * sizeof(double);
            for (RT::List<Channel>::iterator i = channels.begin(), end = channels.end(); i != end; ++i)
                if (i->block)
                    data[n++] = i->block->getValue(i->type, i->index);


	    //Check tooBig as in KZ
	    //cout<<"token"<<token<<endl;
	    if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(sizeof(data))) {
	      fifo.write(&token,sizeof(token));
	      fifo.write(data,sizeof(data));
	    }
	    

	    //DataRecorder::Panel::AppendMe(data);
	    //Commented in original
            //fifo.write(&token, sizeof(token));
            //fifo.write(data, sizeof(data));
        }
    count++;
    counter %= downsample_rate;
}
*/


// Event handler
void DataRecorder::Panel::receiveEvent(const Event::Object *event)
{
    if (event->getName() == Event::IO_BLOCK_INSERT_EVENT)
        {
            IO::Block *block = reinterpret_cast<IO::Block *> (event->getParam("block"));
            blockPtrList.push_back(block);
            blockList->addItem(QString::fromStdString(block->getName()) + " " + QString::number(block->getID()));
            buildChannelList();
        }
    else if (event->getName() == Event::IO_BLOCK_REMOVE_EVENT)
        {
            IO::Block *block = reinterpret_cast<IO::Block *> (event->getParam("block"));
            QString name = QString::fromStdString(block->getName()) + " " + QString::number(block->getID());
            int n = 0;
            for (; n < blockList->count() && blockList->itemText(n) != name; ++n) ;
            if (n < blockList->count()) 
                blockList->removeItem(n);
            blockPtrList.erase(blockPtrList.begin() + n);

            for (RT::List<Channel>::iterator i = channels.begin(), end = channels.end(); i != end; ++i)
                if (i->block == block) {
                    if (recording) i->block = 0;
                    RemoveChannelEvent RTevent(recording, channels, *i);
                    if (!RT::System::getInstance()->postEvent(&RTevent)) {
                        QList<QListWidgetItem*> channelItems = selectionBox->findItems(i->name, Qt::MatchExactly);
                        if (!channelItems.isEmpty()) {
                            /* Use takeItem(row) to remove the channel item. */
                            selectionBox->takeItem(selectionBox->row(channelItems.takeFirst()));
                        } 
                    }
                }
            buildChannelList();
        }
    else if (event->getName() == Event::OPEN_FILE_EVENT)
        {
            QString filename(reinterpret_cast<char*> (event->getParam("filename")));
            OpenFileEvent RTevent(filename, fifo);
            RT::System::getInstance()->postEvent(&RTevent);
        }
    else if (event->getName() == Event::START_RECORDING_EVENT)
        {
            StartRecordingEvent RTevent(recording, fifo);
            RT::System::getInstance()->postEvent(&RTevent);
        }
    else if (event->getName() == Event::STOP_RECORDING_EVENT)
        {
            StopRecordingEvent RTevent(recording, fifo);
            RT::System::getInstance()->postEvent(&RTevent);
        }
    else if (event->getName() == Event::ASYNC_DATA_EVENT)
        {
            AsyncDataEvent RTevent(reinterpret_cast<double *> (event->getParam("data")),*reinterpret_cast<size_t *> (event->getParam("size")), fifo);
            RT::System::getInstance()->postEvent(&RTevent);
        }
    else if( event->getName() == Event::RT_POSTPERIOD_EVENT )
        {
            sleep.tv_nsec = RT::System::getInstance()->getPeriod(); // Update recording thread sleep time
            buildChannelList();
        }
}

// RT Event Handler
void DataRecorder::Panel::receiveEventRT(const Event::Object *event)
{
    if (event->getName() == Event::OPEN_FILE_EVENT)
        {
            QString filename = QString(reinterpret_cast<char*> (event->getParam("filename")));
            data_token_t token;
            token.type = DataRecorder::OPEN;
            token.size = filename.length() + 1;
            token.time = RT::OS::getTime();

	    //Added as in KZ but OPEN_FILE_EVENT was actually missing in KZ
	    if (!fifo.tooBig(sizeof(token))){
	      fifo.write(&token, sizeof(token));
	      fifo.write(filename.toLatin1().constData(), token.size);
	    }

	    //Commented from original
            // fifo.write(&token, sizeof(token));
            // fifo.write(filename.toLatin1().constData(), token.size);
        }
    else if (event->getName() == Event::START_RECORDING_EVENT)
        {
            data_token_t token;
            recording = true;
            token.type = DataRecorder::START;
            token.size = 0;
            token.time = RT::OS::getTime();

	    //Added from KZ 
	    if (!fifo.tooBig(sizeof(token)))
	      fifo.write(&token,sizeof(token));

	    //Commented from original
            //fifo.write(&token, sizeof(token));
        }
    else if (event->getName() == Event::STOP_RECORDING_EVENT)
        {
            data_token_t token;
            recording = false;
            token.type = DataRecorder::STOP;
            token.size = 0;
            token.time = RT::OS::getTime();

	    //Added from KZ 
	    if (!fifo.tooBig(sizeof(token)))
	      fifo.write(&token,sizeof(token));

	    //Commented from original
            //fifo.write(&token, sizeof(token));
        }
    else if (event->getName() == Event::ASYNC_DATA_EVENT)
        {
            size_t size = *reinterpret_cast<size_t *> (event->getParam("size"));
            data_token_t token;
            token.type = DataRecorder::ASYNC;
            token.size = size * sizeof(double);
            token.time = RT::OS::getTime();

	    //Added from KZ 
	    if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	      fifo.write(&token,sizeof(token));
	      fifo.write(event->getParam("data"),token.size);
	    }

	    //Commented from original
            // fifo.write(&token, sizeof(token));
            // fifo.write(event->getParam("data"), token.size);
        }
    else if (event->getName() == Event::WORKSPACE_PARAMETER_CHANGE_EVENT)
        {
            data_token_t token;
            token.type = DataRecorder::PARAM;
            token.size = sizeof(param_change_t);
            token.time = RT::OS::getTime();
            param_change_t data;
            data.id = reinterpret_cast<Settings::Object::ID> (event->getParam("object"));
            data.index = reinterpret_cast<size_t> (event->getParam("index"));
            data.step = file.idx;
            data.value = *reinterpret_cast<double *> (event->getParam("value"));

	    //Added from KZ 
	    if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(sizeof(data))) {
	      fifo.write(&token,sizeof(token));
	      fifo.write(&data,sizeof(data));
	    }

	    //Commented from original
            // fifo.write(&token, sizeof(token));
            // fifo.write(&data, sizeof(data));
        }
    else if( event->getName() == Event::RT_POSTPERIOD_EVENT )
        {
            sleep.tv_nsec = RT::System::getInstance()->getPeriod(); // Update recording thread sleep time
        }
}

// Populate list of blocks and channels
void DataRecorder::Panel::buildChannelList(void)
{
    channelList->clear();
    if (!blockList->count())
        return;

    // Get block
    IO::Block *block = blockPtrList[blockList->currentIndex()];

    // Get type
    IO::flags_t type;
    switch (typeList->currentIndex())
        {
        case 0:
            type = Workspace::INPUT;
            break;
        case 1:
            type = Workspace::OUTPUT;
            break;
        case 2:
            type = Workspace::PARAMETER;
            break;
        case 3:
            type = Workspace::STATE;
            break;
        case 4:
            type = Workspace::EVENT;
            break;
        default:
            ERROR_MSG("DataRecorder::Panel::buildChannelList : invalid type selection\n");
            typeList->setCurrentIndex(0);
            type = Workspace::INPUT;
        }

    for (size_t i = 0; i < block->getCount(type); ++i)
        channelList->addItem(QString::fromStdString(block->getName(type, i)));

    if(channelList->count())
        rButton->setEnabled(true);
    else
        rButton->setEnabled(false);
}

// Slot for changing data file
void DataRecorder::Panel::changeDataFile(void)
{
    QFileDialog fileDialog(this);
    fileDialog.setFileMode(QFileDialog::AnyFile);
    fileDialog.setWindowTitle("Select Data File");

    QSettings userprefs;
    userprefs.setPath(QSettings::NativeFormat, QSettings::SystemScope, "/usr/local/share/rtxi/");
    fileDialog.setDirectory(userprefs.value("/dirs/data", getenv("HOME")).toString());

    QStringList filterList;
    filterList.push_back("HDF5 files (*.h5)");
    filterList.push_back("All files (*.*)");
    fileDialog.setNameFilters(filterList);
    fileDialog.selectNameFilter("HDF5 files (*.h5)");

    QStringList files;
    if(fileDialog.exec())
        files = fileDialog.selectedFiles();

    QString filename;
    if(files.isEmpty() || files[0] == NULL || files[0] == "/" )
        return;
    else
        filename = files[0];

    if (!filename.toLower().endsWith(QString(".h5")))
        filename += ".h5";

    // Write this directory to the user prefs as most recently used
    userprefs.setValue("/dirs/data", fileDialog.directory().path());

    // Post to event queue
    OpenFileEvent RTevent(filename, fifo);
    RT::System::getInstance()->postEvent(&RTevent);
}

// Insert channel to record into list
void DataRecorder::Panel::insertChannel(void)
{
    if (!blockList->count() || !channelList->count())
        return;

    Channel *channel = new Channel();
    channel->block = blockPtrList[blockList->currentIndex()];
    switch (typeList->currentIndex())
        {
        case 0:
            channel->type = Workspace::INPUT;
            break;
        case 1:
            channel->type = Workspace::OUTPUT;
            break;
        case 2:
            channel->type = Workspace::PARAMETER;
            break;
        case 3:
            channel->type = Workspace::STATE;
            break;
        case 4:
            channel->type = Workspace::EVENT;
            break;
        default:
            ERROR_MSG("DataRecorder::Panel::insertChannel : invalid type selection\n");
            typeList->setCurrentIndex(0);
            channel->type = Workspace::INPUT;
        }
    channel->index = channelList->currentIndex();

    channel->name.sprintf("%s %ld : %s", channel->block->getName().c_str(),
                          channel->block->getID(), channel->block->getName(channel->type, channel->index).c_str());

    if(selectionBox->findItems(QString(channel->name), Qt::MatchExactly).isEmpty())
        {
            InsertChannelEvent RTevent(recording, channels, channels.end(), *channel);
            if (!RT::System::getInstance()->postEvent(&RTevent))
                selectionBox->addItem(channel->name);
        }

    if(selectionBox->count())
        {
            lButton->setEnabled(true);
            if(!fileNameEdit->text().isEmpty())
                {
                    startRecordButton->setEnabled(true);
                }
        }
    else
        {
            startRecordButton->setEnabled(false);
            lButton->setEnabled(false);
        }
}

// Remove channel from recorder list
void DataRecorder::Panel::removeChannel(void)
{
    if(!selectionBox->count() || selectionBox->selectedItems().isEmpty())
        return;

    for (RT::List<Channel>::iterator i = channels.begin(), end = channels.end(); i != end; ++i)
        if (i->name == selectionBox->selectedItems().first()->text())
            {
                RemoveChannelEvent RTevent(recording, channels, *i);
                if (!RT::System::getInstance()->postEvent(&RTevent))
                    selectionBox->takeItem(selectionBox->row(selectionBox->selectedItems().first()));
                break;
            }

    if(selectionBox->count())
        {
            startRecordButton->setEnabled(true);
            lButton->setEnabled(true);
        }
    else
        {
            startRecordButton->setEnabled(false);
            lButton->setEnabled(false);
        }
}

// Register new data tag/stamp
void DataRecorder::Panel::addNewTag(void)
{
    std::string newTag(std::to_string(RT::OS::getTime()));
    newTag += ",";
    newTag += timeStampEdit->text().toStdString();
    dataTags.push_back(newTag);
    timeStampEdit->clear();
    recordStatus->setText("Tagged");
}

// Start recording slot
void DataRecorder::Panel::startRecordClicked(void)
{
    if(fileNameEdit->text().isEmpty())
        {
            QMessageBox::critical(
                this, "Data file not specified.",
                "Please specify a file to write data to.",
                QMessageBox::Ok, QMessageBox::NoButton);
            return;
        }

    StartRecordingEvent RTevent(recording, fifo);
    RT::System::getInstance()->postEvent(&RTevent);
}

// Stop recording slot
void DataRecorder::Panel::stopRecordClicked(void)
{
    StopRecordingEvent RTevent(recording, fifo);
    RT::System::getInstance()->postEvent(&RTevent);
}

// Update downsample rate
void DataRecorder::Panel::updateDownsampleRate(int r)
{
    downsample_rate = r;
}

// Custom event handler
void DataRecorder::Panel::customEvent(QEvent *e)
{
    if (e->type() == QFileExistsEvent)
        {
            mutex.lock();
            CustomEvent * event = static_cast<CustomEvent *>(e);
            FileExistsEventData *data = reinterpret_cast<FileExistsEventData *> (event->getData());
            data->response = QMessageBox::question(this, "File exists",
                                                   "The file already exists. What would you like to do?",
                                                   "Append", "Overwrite", "Cancel", 0, 2);
            recordStatus->setText("Not Recording");
            data->done.wakeAll();
            mutex.unlock();
        }
    else if (e->type() == QSetFileNameEditEvent)
        {
            mutex.lock();
            CustomEvent * event = static_cast<CustomEvent *>(e);
            SetFileNameEditEventData *data = reinterpret_cast<SetFileNameEditEventData *> (event->getData());
            fileNameEdit->setText(data->filename);
            recordStatus->setText("Ready.");
            if(selectionBox->count())
                {
                    startRecordButton->setEnabled(true);
                }
            data->done.wakeAll();
            mutex.unlock();
        }
    else if (e->type() == QDisableGroupsEvent)
        {
            startRecordButton->setEnabled(false);
            stopRecordButton->setEnabled(true);
            closeButton->setEnabled(false);
            channelGroup->setEnabled(false);
            sampleGroup->setEnabled(false);
            recordStatus->setText("Recording...");
        }
    else if (e->type() == QEnableGroupsEvent)
        {
            startRecordButton->setEnabled(true);
            stopRecordButton->setEnabled(false);
            closeButton->setEnabled(true);
            channelGroup->setEnabled(true);
            sampleGroup->setEnabled(true);
            recordStatus->setText("Ready.");
            fileSize->setNum(int(QFile(fileNameEdit->text()).size())/1024.0/1024.0);
            trialLength->setNum(double(RT::System::getInstance()->getPeriod()*1e-9* fixedcount));
            count = 0;
        }
}

void DataRecorder::Panel::doDeferred(const Settings::Object::State &s)
{
    for (int i = 0; i < s.loadInteger("Num Channels"); ++i)
        {
            Channel *channel;
            IO::Block *block;
            std::ostringstream str;
            str << i;

            block	= dynamic_cast<IO::Block *> (Settings::Manager::getInstance()->getObject(s.loadInteger(str.str() + " ID")));
            if (!block)
                continue;

            channel = new Channel();
            channel->block = block;
            channel->type = s.loadInteger(str.str() + " type");
            channel->index = s.loadInteger(str.str() + " index");
            channel->name.sprintf("%s %ld : %s", channel->block->getName().c_str(),
                                  channel->block->getID(), channel->block->getName(channel->type,	channel->index).c_str());

            channels.insert(channels.end(), *channel);
            selectionBox->addItem(channel->name);
            if(selectionBox->count())
                lButton->setEnabled(true);
        }
}

void DataRecorder::Panel::doLoad(const Settings::Object::State &s)
{
    if (s.loadInteger("Maximized"))
        showMaximized();
    else if (s.loadInteger("Minimized"))
        showMinimized();

    downsampleSpin->setValue(s.loadInteger("Downsample"));
    parentWidget()->move(s.loadInteger("X"), s.loadInteger("Y"));
}

void DataRecorder::Panel::doSave(Settings::Object::State &s) const
{
    if (isMaximized())
        s.saveInteger("Maximized", 1);
    else if (isMinimized())
        s.saveInteger("Minimized", 1);

    QPoint pos = parentWidget()->pos();
    s.saveInteger("X", pos.x());
    s.saveInteger("Y", pos.y());

    s.saveInteger("Downsample", downsampleSpin->value());
    s.saveInteger("Num Channels", channels.size());
    size_t n = 0;
    for (RT::List<Channel>::const_iterator i = channels.begin(), end = channels.end(); i != end; ++i)
        {
            std::ostringstream str;
            str << n++;

            s.saveInteger(str.str() + " ID", i->block->getID());
            s.saveInteger(str.str() + " type", i->type);
            s.saveInteger(str.str() + " index", i->index);
        }
}

void *DataRecorder::Panel::bounce(void *param)
{
    Panel *that = reinterpret_cast<Panel *> (param);
    if (that)
        {
            that->processData();
        }
    return 0;
}

void DataRecorder::Panel::processData(void)
{
    enum
    {
        CLOSED, OPENED, RECORD,
    } state = CLOSED;

    tokenRetrieved = false;
    
    for (;;)
        {
	  if(!tokenRetrieved)
	    {
	      // Returns true if data was available and retrieved
	      if(fifo.read(&_token, sizeof(_token)))
		tokenRetrieved = true;
	      else
		{
		  // Sleep loop then restart if no token was retrieved
		  nanosleep(&sleep, NULL);
		  continue;
		}
	    }


	    //**********************************************************************************
	    //*                               Original SYNC                                    *
	    //**********************************************************************************
            // if (_token.type == SYNC)
	    //   {
		
	    // 	if (state == RECORD)
	    // 	  {
	    // 	    double data[_token.size / sizeof(double)];
	    // 	    if(!fifo.read(data, _token.size))
	    // 	      continue; // Restart loop if data is not available
	    // 	    H5PTappend(file.cdata, 1, data);
	    // 	    ++file.idx;
	    // 	  }

	    //   }


	  
            if (_token.type == _START)
	      {
		if (state == RECORD)
		  {
		    //Only for continuous recording like in SAP
		    time(&current_time);
		    ss << current_time; 
		    file_id = ss.str();
		    file_name = file_path + "file_" + file_id + ".txt";
		    //file_name_flt = file_path_flt + "file_" + file_id + "_flt.txt";
		    file_to_write.open(file_name);
		    //file_to_write_flt.open(file_name_flt);
		    if(!file_to_write.is_open()){
		      cout<<"Unable to open file in _START/RECORD. File is open:  "<<file_to_write.is_open()<<endl;
		      if (getcwd(cwd, sizeof(cwd)) != NULL) {
			cout<<"Current working dir: "<<cwd<<endl;
		      } 
		    }
		    cout<<"File opened"<<endl;

		    // if(!file_to_write_flt.is_open()){
		    //   cout<<"Unable to open file in _START/RECORD. File is open:  "<<file_to_write_flt.is_open()<<endl;
		    //   if (getcwd(cwd, sizeof(cwd)) != NULL) {
		    // 	cout<<"Current working dir: "<<cwd<<endl;
		    //   } 
		    // }
		    // cout<<"File flt opened"<<endl;

		    
		  }
	      }

            else if (_token.type == _STOP)
	      {
		if (state == RECORD)
		  {
		    file_to_write.close();
		    cout<<"File written"<<endl;
		    //file_to_write_flt.close();
		    //cout<<"File flt written"<<endl;
		    ss.str(string());
  
		  }
	      }

            else if (_token.type == SYNC)
	      {
		if (state == RECORD)
		  {
		    double data[_token.size / sizeof(double)];
		    if(!fifo.read(data, _token.size))
		      {
			continue; // Restart loop if data is not available
		      }
		    //End original
		    if(file_to_write.is_open()){
		      file_to_write << data[0]<<'\n'; //data
		      //file_to_write_flt << data[1]<<'\n'; //data_flt
		    }
		    else
		      {
			cout<<"File or file_flt is not opened in SYNC/RECORD. File is open:  "<<file_to_write.is_open()<<endl;
			if (getcwd(cwd, sizeof(cwd)) != NULL) {
			  cout<<"Current working dir: "<<cwd<<endl;
			}
		      }

		    //Append once in a second data to H5 file
		    if(h5_counter==max_h5_counter)
		      {
			H5PTappend(file.cdata, 1, data);
			++file.idx;
			h5_counter=0;
		      }
		    h5_counter++;

		  }
	      }
            else if (_token.type == ASYNC)
                {
                    if (state == RECORD)
                        {
                            double data[_token.size / sizeof(double)];
                            if(!fifo.read(data, _token.size))
                                continue; // Restart loop if data is not available
                            if (data)
                                {
                                    hsize_t array_size[] = { _token.size / sizeof(double) };
                                    hid_t array_space = H5Screate_simple(1, array_size,	array_size);
                                    hid_t array_type = H5Tarray_create(H5T_IEEE_F64LE, 1,	array_size);

                                    QString data_name = QString::number(static_cast<unsigned long long> (_token.time));
                                    hid_t adata = H5Dcreate(file.adata, data_name.toLatin1().constData(),
                                                            array_type, array_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                                    H5Dwrite(adata, array_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

                                    H5Dclose(adata);
                                    H5Tclose(array_type);
                                    H5Sclose(array_space);
                                }
                        }
                }
            else if (_token.type == OPEN)
                {
                    if (state == RECORD)
                        stopRecording(_token.time);
                    if (state != CLOSED)
                        closeFile();
                    char filename_string[_token.size];
                    if(!fifo.read(filename_string, _token.size))
                        continue; // Restart loop if data is not available
                    QString filename = filename_string;
		    //cout<<"Going to open file"<<endl;
                    if (openFile(filename))
                        state = CLOSED;
                    else
                        state = OPENED;
                }
            else if (_token.type == CLOSE)
                {
                    if (state == RECORD)
                        stopRecording(RT::OS::getTime());
                    if (state != CLOSED)
                        closeFile();
                    state = CLOSED;
                }
            else if (_token.type == START)
                {
                    if (state == OPENED)
                        {
                            count = 0;
                            startRecording(_token.time);
                            state = RECORD;
                            QEvent *event = new QEvent(static_cast<QEvent::Type>QDisableGroupsEvent);
                            QApplication::postEvent(this, event);
                        }
                }
            else if (_token.type == STOP)
                {
                    if (state == RECORD)
                        {
                            stopRecording(_token.time);
                            state = OPENED;
                            fixedcount = count;
                            QEvent *event = new QEvent(static_cast<QEvent::Type>QEnableGroupsEvent);
                            QApplication::postEvent(this, event);
                        }
                }
            else if (_token.type == DONE)
                {
                    if (state == RECORD)
                        stopRecording(_token.time);
                    if (state != CLOSED)
                        closeFile(true);
                    break;
                }
            else if (_token.type == PARAM)
                {
                    param_change_t data;
                    if(!fifo.read(&data, sizeof(data)))
                        continue; // Restart loop if data is not available

                    IO::Block	*block = dynamic_cast<IO::Block *> (Settings::Manager::getInstance()->getObject(data.id));

                    if (block && state == RECORD)
                        {
                            param_hdf_t param = { data.step, data.value, };

                            hid_t param_type;
                            param_type = H5Tcreate(H5T_COMPOUND, sizeof(param_hdf_t));
                            H5Tinsert(param_type, "index", HOFFSET(param_hdf_t,index), H5T_STD_I64LE);
                            H5Tinsert(param_type, "value", HOFFSET(param_hdf_t,value), H5T_IEEE_F64LE);

                            QString parameter_name = QString::number(block->getID()) + " "
                                                     + QString::fromStdString(block->getName()) + " : "
                                                     + QString::fromStdString(block->getName(Workspace::PARAMETER, data.index));

                            hid_t data = H5PTopen(file.pdata, parameter_name.toLatin1().constData());
                            H5PTappend(data, 1, &param);
                            H5PTclose(data);
                            H5Tclose(param_type);
                        }
                }
            tokenRetrieved = false;
        }
}

int DataRecorder::Panel::openFile(QString &filename)
{
#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread))
        {
            ERROR_MSG("DataRecorder::Panel::openFile : called by invalid thread\n");
            PRINT_BACKTRACE();
        }
#endif

    if (QFile::exists(filename))
        {
            mutex.lock();
            CustomEvent *event = new CustomEvent(static_cast<QEvent::Type>QFileExistsEvent);
            FileExistsEventData data;
            event->setData(static_cast<void *>(&data));
            data.filename = filename;
            QApplication::postEvent(this, event);
            data.done.wait(&mutex);
            mutex.unlock();

            if (data.response == 0)   // append
                {
                    file.id = H5Fopen(filename.toLatin1().constData(), H5F_ACC_RDWR, H5P_DEFAULT);
                    size_t trial_num;
                    QString trial_name;
                    H5Eset_auto(H5E_DEFAULT, NULL, NULL);
                    for (trial_num = 1;; ++trial_num)
                        {
                            trial_name = "/Trial" + QString::number(trial_num);
                            file.trial = H5Gopen(file.id, trial_name.toLatin1().constData(), H5P_DEFAULT);
                            if (file.trial < 0)
                                {
                                    H5Eclear(H5E_DEFAULT);
                                    break;
                                }
                            else
                                {
                                    H5Gclose(file.trial);
                                }
                        }
                    trialNum->setNum(int(trial_num)-1);
                }
            else if (data.response == 1)     //overwrite
                {
                    file.id = H5Fcreate(filename.toLatin1().constData(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                    trialNum->setText("0");
                }
            else
                {
                    return -1;
                }
        }
    else
        {
            file.id = H5Fcreate(filename.toLatin1().constData(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            trialNum->setText("0");
        }
    if (file.id < 0)
        {
            H5E_type_t error_type;
            size_t error_size;
            error_size = H5Eget_msg(file.id, &error_type, NULL, 0);
            char error_msg[error_size + 1];
            H5Eget_msg(file.id, &error_type, error_msg, error_size);
            error_msg[error_size] = 0;
            H5Eclear(file.id);

            ERROR_MSG("DataRecorder::Panel::processData : failed to open \"%s\" for writing with error : %s\n", filename.toStdString().c_str(),error_msg);
            return -1;
        }

    mutex.lock();
    CustomEvent *event = new CustomEvent(static_cast<QEvent::Type>QSetFileNameEditEvent);
    SetFileNameEditEventData data;
    data.filename = filename;
    event->setData(static_cast<void*>(&data));
    QApplication::postEvent(this, event);
    data.done.wait(&mutex);
    mutex.unlock();

    return 0;
}

void DataRecorder::Panel::closeFile(bool shutdown)
{
#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread))
        {
            ERROR_MSG("DataRecorder::Panel::closeFile : called by invalid thread\n");
            PRINT_BACKTRACE();
        }
#endif

				if(!dataTags.empty())
				{
    // Write tags to data file
    hid_t tag_type, tag_space, data;
    herr_t status;
    hsize_t dims[1] = {1};
    tag_type = H5Tcreate(H5T_STRING, TAG_SIZE);
    tag_space = H5Screate_simple(1, dims, NULL);

    // Create group for tags
    file.tdata = H5Gcreate(file.id, "Tags", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Iterate over vector (buffer) and put into data file
    size_t i = 0;
    for(std::vector<std::string>::iterator it = dataTags.begin(); it != dataTags.end(); ++it)
        {
            data = H5Dcreate(file.tdata, std::string("Tag " + std::to_string(i++)).c_str(), tag_type, tag_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(data, tag_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, it->c_str());
        }
    dataTags.clear();

    // Close all open structs
    H5Dclose(data);
    H5Sclose(tag_space);
    H5Tclose(tag_type);
    H5Gclose(file.tdata);
				}

    // Close file
    H5Fclose(file.id);

    if (!shutdown)
        {
            mutex.lock();
            CustomEvent *event = new CustomEvent(static_cast<QEvent::Type>QSetFileNameEditEvent);
            SetFileNameEditEventData data;
            data.filename = "";
            event->setData(static_cast<void*>(&data));
            QApplication::postEvent(this, event);
            data.done.wait(&mutex);
            mutex.unlock();
        }
}


int DataRecorder::Panel::startRecording(long long timestamp)
{
#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread))
        {
            ERROR_MSG("DataRecorder::Panel::startRecording : called by invalid thread\n");
            PRINT_BACKTRACE();
        }
#endif

    size_t trial_num;
    QString trial_name;

    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    for (trial_num = 1;; ++trial_num)
        {
            trial_name = "/Trial" + QString::number(trial_num);
            file.trial = H5Gopen(file.id, trial_name.toLatin1().constData(), H5P_DEFAULT);

            if (file.trial < 0)
                {
                    H5Eclear(H5E_DEFAULT);
                    break;
                }
            else
                H5Gclose(file.trial);
        }

    trialNum->setNum(int(trial_num));
    file.trial = H5Gcreate(file.id, trial_name.toLatin1().constData(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    file.pdata = H5Gcreate(file.trial, "Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    file.adata = H5Gcreate(file.trial, "Asynchronous Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    file.sdata = H5Gcreate(file.trial, "Synchronous Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    file.sysdata = H5Gcreate(file.trial, "System Settings", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t string_type = H5Tcopy(H5T_C_S1);
    size_t string_size = 1024;
    H5Tset_size(string_type, string_size);
    hid_t data;

    data = H5Dcreate(file.trial, "Version", string_type,
                     scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		std::string version_string = QString(VERSION).toStdString();
		char * version_c_string = new char[version_string.length()+1];
		std::strcpy(version_c_string, version_string.c_str());
    H5Dwrite(data, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, version_c_string);
    delete[] version_c_string;
    H5Dclose(data);

    long long period = RT::System::getInstance()->getPeriod();
    data = H5Dcreate(file.trial, "Period (ns)", H5T_STD_U64LE, scalar_space,
                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(data, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &period);
    H5Dclose(data);

    long long downsample = downsample_rate;
    data = H5Dcreate(file.trial, "Downsampling Rate", H5T_STD_U64LE,
                     scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(data, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &downsample);
    H5Dclose(data);

    data = H5Dcreate(file.trial, "Timestamp Start (ns)", H5T_STD_U64LE,
                     scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(data, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &timestamp);
    H5Dclose(data);

    data = H5Dcreate(file.trial, "Date", string_type,
                     scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		std::string date_string = QDateTime::currentDateTime().toString(Qt::ISODate).toStdString();
		char * date_c_string = new char[date_string.length()+1];
		std::strcpy(date_c_string, date_string.c_str());
    H5Dwrite(data, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, date_c_string);
    delete[] date_c_string;
    H5Dclose(data);

    hid_t param_type;
    param_type = H5Tcreate(H5T_COMPOUND, sizeof(param_hdf_t));
    H5Tinsert(param_type, "index", HOFFSET(param_hdf_t,index), H5T_STD_I64LE);
    H5Tinsert(param_type, "value", HOFFSET(param_hdf_t,value), H5T_IEEE_F64LE);

    for (RT::List<Channel>::iterator i = channels.begin(), end = channels.end(); i != end; ++i)
        {
            IO::Block *block = i->block;
            for (size_t j = 0; j < block->getCount(Workspace::PARAMETER); ++j)
                {
                    QString parameter_name = QString::number(block->getID()) + " "
                                             + QString::fromStdString(block->getName()) + " : " + QString::fromStdString(block->getName(Workspace::PARAMETER, j));
                    data = H5PTcreate_fl(file.pdata, parameter_name.toLatin1().constData(),	param_type, sizeof(param_hdf_t), -1);
                    struct param_hdf_t value = { 0, block->getValue(Workspace::PARAMETER, j),};
                    H5PTappend(data, 1, &value);
                    H5PTclose(data);
                }
            for (size_t j = 0; j < block->getCount(Workspace::COMMENT); ++j)
                {
                    QString comment_name = QString::number(block->getID()) + " "
                                           + QString::fromStdString(block->getName()) + " : " + QString::fromStdString(block->getName(Workspace::COMMENT, j));
                    hsize_t	dims = dynamic_cast<Workspace::Instance *> (block)->getValueString(Workspace::COMMENT, j).size() + 1;
                    hid_t comment_space = H5Screate_simple(1, &dims, &dims);
                    data = H5Dcreate(file.pdata, comment_name.toLatin1().constData(), H5T_C_S1,	comment_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    H5Dwrite(data, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT,	dynamic_cast<Workspace::Instance *> (block)->getValueString(Workspace::COMMENT, j).c_str());
                    H5Dclose(data);
                }
        }

    H5Tclose(param_type);

    //Added :: to ispunct at line 2094 below so that ispunct is taken from the main namespace instead of std namespace.
    size_t count = 0;
    for (RT::List<Channel>::iterator i = channels.begin(), end = channels.end(); i != end; ++i)
        {
            std::string rec_chan_name = std::to_string(++count) + " " + i->name.toStdString();
            rec_chan_name.erase(std::remove_if(rec_chan_name.begin(), rec_chan_name.end(), &::ispunct), rec_chan_name.end());
            hid_t data = H5Dcreate(file.sdata, rec_chan_name.c_str(), string_type, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            H5Dwrite(data, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, rec_chan_name.c_str());
            H5Dclose(data);
        }

    DAQ::Device *dev;
    {
        struct find_daq_t info = { 0, 0, };
        DAQ::Manager::getInstance()->foreachDevice(findDAQDevice, &info);
        dev = info.device;
    }

				// Save channel configurations
				if(dev)
					for(size_t i=0; i<dev->getChannelCount(DAQ::AI); ++i)
						if(dev->getChannelActive(DAQ::AI,static_cast<DAQ::index_t>(i)))
						{
							std::string chan_name = "Analog Channel " + std::to_string(i);
							file.chandata = H5Gcreate(file.sysdata, chan_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

							hid_t data = H5Dcreate(file.chandata, "Range", string_type, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
							std::string range_string = dev->getAnalogRangeString(DAQ::AI,static_cast<DAQ::index_t>(i),dev->getAnalogRange(DAQ::AI,static_cast<DAQ::index_t>(i)));
							char * range_c_string = new char[range_string.length()+1];
							std::strcpy(range_c_string, range_string.c_str());
							H5Dwrite(data, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, range_c_string);
							delete[] range_c_string;

							data = H5Dcreate(file.chandata, "Reference", string_type, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
							std::string ref_string = dev->getAnalogReferenceString(DAQ::AI,static_cast<DAQ::index_t>(i),dev->getAnalogReference(DAQ::AI,static_cast<DAQ::index_t>(i)));
							char * ref_c_string = new char[ref_string.length()+1];
							std::strcpy(ref_c_string, ref_string.c_str());
							H5Dwrite(data, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, ref_c_string);
							delete[] ref_c_string;

							double scale = dev->getAnalogGain(DAQ::AI,static_cast<DAQ::index_t>(i));
							data = H5Dcreate(file.chandata, "Gain", H5T_IEEE_F64LE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
							H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &scale); 

							double offset = dev->getAnalogZeroOffset(DAQ::AI,static_cast<DAQ::index_t>(i));
							data = H5Dcreate(file.chandata, "Offset", H5T_IEEE_F64LE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
							H5Dwrite(data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &offset); 

							int downsample = dev->getAnalogDownsample(DAQ::AI,static_cast<DAQ::index_t>(i));
							data = H5Dcreate(file.chandata, "Downsample", H5T_STD_I16LE, scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
							H5Dwrite(data, H5T_STD_I16LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &downsample);
							H5Dclose(data);
						}

    H5Tclose(string_type);
    H5Sclose(scalar_space);

    if (channels.size())
        {
            hsize_t array_size[] = { channels.size() };
            hid_t array_type = H5Tarray_create(H5T_IEEE_F64LE, 1, array_size);
            file.cdata = H5PTcreate_fl(file.sdata, "Channel Data", array_type, (hsize_t) 64, 1);
            H5Tclose(array_type);
        }

    file.idx = 0;

    return 0;
}


void DataRecorder::Panel::stopRecording(long long timestamp)
{
#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread))
        {
            ERROR_MSG("DataRecorder::Panel::stopRecording : called by invalid thread\n");
            PRINT_BACKTRACE();
        }
#endif

    // Write stop time to data file
    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t data = H5Dcreate(file.trial, "Timestamp Stop (ns)", H5T_STD_U64LE,
                           scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(data, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &timestamp);
    H5Dclose(data);

    // Write trial length to data file
    fixedcount = count;
    long long period = RT::System::getInstance()->getPeriod();
    long long datalength = period * fixedcount;
    data = H5Dcreate(file.trial, "Trial Length (ns)", H5T_STD_U64LE,
                     scalar_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(data, H5T_STD_U64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datalength);
    H5Dclose(data);

    // Close all open structs
    H5Sclose(scalar_space);
    H5PTclose(file.cdata);
    H5Gclose(file.sdata);
    H5Gclose(file.pdata);
    H5Gclose(file.adata);
    H5Gclose(file.sysdata);
    H5Gclose(file.chandata);
    H5Gclose(file.trial);

    H5Fflush(file.id, H5F_SCOPE_LOCAL);
    void *file_handle;
    H5Fget_vfd_handle(file.id, H5P_DEFAULT, &file_handle);
    if (fsync(*static_cast<int *> (file_handle)))
        {
            DEBUG_MSG("DataRecorder::Panel::stopRecording : fsync failed, running sync\n");
            sync();
        }
}

extern "C" Plugin::Object *createRTXIPlugin(void *)
{
    return DataRecorder::Plugin::getInstance();
}

DataRecorder::Plugin::Plugin(void)
{
    // get the HDF data recorder buffer size from user preference
    QSettings userprefs;
    userprefs.setPath(QSettings::NativeFormat, QSettings::SystemScope, "/usr/local/share/rtxi/");
    buffersize = (userprefs.value("/system/HDFbuffer", 10).toInt())*20971520; //1048576
    MainWindow::getInstance()->createSystemMenuItem("Data Recorder",this,SLOT(createDataRecorderPanel(void)));
}

DataRecorder::Plugin::~Plugin(void)
{
    while (panelList.size())
        delete panelList.front();
    instance = 0;
}

DataRecorder::Panel *DataRecorder::Plugin::createDataRecorderPanel(void)
{
    Panel *panel = new Panel(MainWindow::getInstance()->centralWidget(), buffersize);
    panelList.push_back(panel);
    return panel;
}

void DataRecorder::Plugin::removeDataRecorderPanel(DataRecorder::Panel *panel)
{
    panelList.remove(panel);
}

void DataRecorder::Plugin::doDeferred(const Settings::Object::State &s)
{
    size_t i = 0;
    for (std::list<Panel *>::iterator j = panelList.begin(), end = panelList.end(); j != end; ++j)
        (*j)->deferred(s.loadState(QString::number(i++).toStdString()));
}

void DataRecorder::Plugin::doLoad(const Settings::Object::State &s)
{
    for (size_t i = 0; i < static_cast<size_t> (s.loadInteger("Num Panels")); ++i)
        {
            Panel *panel = new Panel(MainWindow::getInstance()->centralWidget(), buffersize);
            panelList.push_back(panel);
            panel->load(s.loadState(QString::number(i).toStdString()));
        }
}

void DataRecorder::Plugin::doSave(Settings::Object::State &s) const
{
    s.saveInteger("Num Panels", panelList.size());
    size_t n = 0;
    for (std::list<Panel *>::const_iterator i = panelList.begin(), end = panelList.end(); i != end; ++i)
        s.saveState(QString::number(n++).toStdString(), (*i)->save());
}

static Mutex mutex;
DataRecorder::Plugin *DataRecorder::Plugin::instance = 0;

DataRecorder::Plugin *DataRecorder::Plugin::getInstance(void)
{
    if (instance)
        return instance;

    /*************************************************************************
     * Seems like alot of hoops to jump through, but allocation isn't        *
     *   thread-safe. So effort must be taken to ensure mutual exclusion.    *
     *************************************************************************/

    Mutex::Locker lock(&::mutex);
    if (!instance)
        instance = new Plugin();

    return instance;
}
