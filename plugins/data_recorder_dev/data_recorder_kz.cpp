/*
 * Copyright (C) 2004 Boston University
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License as
 *  published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


// Code for Spectral Analysis begins on line 488.

#include <qapplication.h> 
#include <qcombobox.h>
#include <qfiledialog.h>
#include <qgroupbox.h>
#include <qhbox.h>
#include <qhgroupbox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qlistbox.h>
#include <qmessagebox.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qstringlist.h>
#include <qvbox.h>
#include <qvgroupbox.h>
#include <qwaitcondition.h>

#include <compiler.h>
#include <debug.h>
#include <main_window.h>
#include <sstream>
#include <workspace.h>
#include <data_recorder.h>

#include <qpixmap.h>
#include <qpainter.h>
#include <qtimer.h>
#include <qvalidator.h>
#include <cmath>
#include <main_window.h>
#include <rt.h>

#define QFileExistsEvent            (QEvent::User+0)
#define QNoFileOpenEvent            (QEvent::User+1)
#define QSetFileNameEditEvent       (QEvent::User+2)
#define QDisableGroupsEvent         (QEvent::User+3)
#define QEnableGroupsEvent          (QEvent::User+4)

struct param_hdf_t {
    long long index;
    double value;
};

namespace {

    void buildBlockPtrList(IO::Block *block,void *arg) {
        std::vector<IO::Block *> *list = reinterpret_cast<std::vector<IO::Block *> *>(arg);
        list->push_back(block);
    };

    struct FileExistsEventData {
        QString filename;
        int response;
        QWaitCondition done;
    };

    struct SetFileNameEditEventData {
        QString filename;
        QWaitCondition done;
    };

    class InsertChannelEvent : public RT::Event {

    public:

        InsertChannelEvent(bool &,RT::List<DataRecorder::Channel> &,RT::List<DataRecorder::Channel>::iterator,DataRecorder::Channel &);
        ~InsertChannelEvent(void);

        int callback(void);

    private:

        bool &recording;
        RT::List<DataRecorder::Channel> &channels;
        RT::List<DataRecorder::Channel>::iterator end;
        DataRecorder::Channel &channel;

    }; // class InsertChannelEvent

    class RemoveChannelEvent : public RT::Event {

    public:

        RemoveChannelEvent(bool &,RT::List<DataRecorder::Channel> &,DataRecorder::Channel &);
        ~RemoveChannelEvent(void);

        int callback(void);

    private:

        bool &recording;
        RT::List<DataRecorder::Channel> &channels;
        DataRecorder::Channel &channel;

    }; // class RemoveChannelEvent

    class OpenFileEvent : public RT::Event {

    public:

        OpenFileEvent(QString &,Fifo &);
        ~OpenFileEvent(void);

        int callback(void);

    private:

        QString &filename;
        Fifo &fifo;

    }; // class OpenFileEvent

    class StartRecordingEvent : public RT::Event {

    public:

        StartRecordingEvent(bool &,Fifo &);
        ~StartRecordingEvent(void);

        int callback(void);

    private:

        bool &recording;
        Fifo &fifo;

    }; // class StartRecordingEvent

    class StopRecordingEvent : public RT::Event {

    public:

        StopRecordingEvent(bool &,Fifo &);
        ~StopRecordingEvent(void);

        int callback(void);

    private:

        bool &recording;
        Fifo &fifo;

    }; //class StopRecordingEvent

    class AsyncDataEvent : public RT::Event {

    public:

        AsyncDataEvent(const double *,size_t,Fifo &);
        ~AsyncDataEvent(void);

        int callback(void);

    private:

        const double *data;
        size_t size;
        Fifo &fifo;

    }; // class AsyncDataEvent

    class DoneEvent : public RT::Event {

    public:

        DoneEvent(Fifo &);
        ~DoneEvent(void);

        int callback(void);

    private:

        Fifo &fifo;

    }; // class DoneEvent

}; // namespace

InsertChannelEvent::InsertChannelEvent(bool &r,RT::List<DataRecorder::Channel> & l,RT::List<DataRecorder::Channel>::iterator e,DataRecorder::Channel &c)
    : recording(r), channels(l), end(e), channel(c) {}

InsertChannelEvent::~InsertChannelEvent(void) {}

int InsertChannelEvent::callback(void) {
    if(recording)
        return -1;

    channels.insertRT(end,channel);
    return 0;
}

RemoveChannelEvent::RemoveChannelEvent(bool &r,RT::List<DataRecorder::Channel> & l,DataRecorder::Channel &c)
    : recording(r), channels(l), channel(c) {}

RemoveChannelEvent::~RemoveChannelEvent(void) {}

int RemoveChannelEvent::callback(void) {
    if(recording)
        return -1;

    channels.removeRT(channel);
    return 0;
}

OpenFileEvent::OpenFileEvent(QString &n,Fifo &f)
    : filename(n), fifo(f) {}

OpenFileEvent::~OpenFileEvent(void) {}

int OpenFileEvent::callback(void) {
    DataRecorder::data_token_t token;

    token.type = DataRecorder::OPEN;
    token.size = filename.length()+1;
    token.time = RT::OS::getTime();

    if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	    fifo.write(&token,sizeof(token));
	    fifo.write(filename.latin1(),token.size);
    }

    return 0;
}

StartRecordingEvent::StartRecordingEvent(bool &r,Fifo &f)
    : recording(r), fifo(f) {}

StartRecordingEvent::~StartRecordingEvent(void) {}

int StartRecordingEvent::callback(void) {
    DataRecorder::data_token_t token;

    recording = true;

    token.type = DataRecorder::START;
    token.size = 0;
    token.time = RT::OS::getTime();

    if (!fifo.tooBig(sizeof(token)))
	    fifo.write(&token,sizeof(token));

    return 0;
}

StopRecordingEvent::StopRecordingEvent(bool &r,Fifo &f)
    : recording(r), fifo(f) {}

StopRecordingEvent::~StopRecordingEvent(void) {}

int StopRecordingEvent::callback(void) {
    DataRecorder::data_token_t token;

    recording = false;

    token.type = DataRecorder::STOP;
    token.size = 0;
    token.time = RT::OS::getTime();

	if (!fifo.tooBig(sizeof(token)))
	    fifo.write(&token,sizeof(token));

    return 0;
}

AsyncDataEvent::AsyncDataEvent(const double *d,size_t s,Fifo &f)
    : data(d), size(s), fifo(f) {}

AsyncDataEvent::~AsyncDataEvent(void) {}

int AsyncDataEvent::callback(void) {
    DataRecorder::data_token_t token;

    token.type = DataRecorder::ASYNC;
    token.size = size*sizeof(double);
    token.time = RT::OS::getTime();

	if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	    fifo.write(&token,sizeof(token));
	    fifo.write(data,token.size);
	}

    return 1;
}

DoneEvent::DoneEvent(Fifo &f)
    : fifo(f) {}

DoneEvent::~DoneEvent(void) {}

int DoneEvent::callback(void) {
    DataRecorder::data_token_t token;

    token.type = DataRecorder::DONE;
    token.size = 0;
    token.time = RT::OS::getTime();

	if (!fifo.tooBig(sizeof(token)))
	    fifo.write(&token,sizeof(token));

    return 0;
}

void DataRecorder::startRecording(void) {
    Event::Object event(Event::START_RECORDING_EVENT);

    if(RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

void DataRecorder::stopRecording(void) {
    Event::Object event(Event::STOP_RECORDING_EVENT);

    if(RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

void DataRecorder::postAsyncData(const double *data,size_t size) {

    Event::Object event(Event::ASYNC_DATA_EVENT);
    event.setParam("data",const_cast<double *>(data));
    event.setParam("size",&size);

    if(RT::OS::isRealtime())
        Event::Manager::getInstance()->postEventRT(&event);
    else
        Event::Manager::getInstance()->postEvent(&event);
}

DataRecorder::Channel::Channel(void) {}

DataRecorder::Channel::~Channel(void) {}

static IO::channel_t chans[] = {
    {
        "Me",
        "",
        IO::OUTPUT,
    },
};
static size_t num_chans = sizeof(chans)/sizeof(chans[0]);

DataRecorder::Panel::Panel(QWidget *parent)
    : QWidget(parent,0,Qt::WStyle_NormalBorder|Qt::WDestructiveClose),
      RT::Thread(RT::Thread::MinimumPriority), fifo(1048576), recording(false), IO::Block("Connect",chans,num_chans) {
    setCaption(QString::number(getID())+" Data Recorder");

    QHBox *hbox;
    QVBox *vbox;

    QBoxLayout *layout = new QVBoxLayout(this);
    layout->setMargin(5);

    channelBox = new QHGroupBox("Channel Selection",this);
    channelBox->setInsideMargin(5);
    layout->addWidget(channelBox);

    vbox = new QVBox(channelBox);
    vbox->setMaximumHeight(125);

    hbox = new QHBox(vbox);
    (new QLabel("Block:",hbox))->setFixedWidth(75);
    blockList = new QComboBox(hbox);
    blockList->setFixedWidth(150);
    QObject::connect(blockList,SIGNAL(activated(int)),this,SLOT(buildChannelList(void)));

    hbox = new QHBox(vbox);
    (new QLabel("Type:",hbox))->setFixedWidth(75);
    typeList = new QComboBox(hbox);
    typeList->setFixedWidth(150);
    typeList->insertItem("Input");
    typeList->insertItem("Output");
    typeList->insertItem("Parameter");
    typeList->insertItem("State");
    typeList->insertItem("Event");
    QObject::connect(typeList,SIGNAL(activated(int)),this,SLOT(buildChannelList(void)));

    hbox = new QHBox(vbox);
    (new QLabel("Channel:",hbox))->setFixedWidth(75);
    channelList = new QComboBox(hbox);
    channelList->setFixedWidth(150);

    vbox = new QVBox(channelBox);
    vbox->setMaximumHeight(100);

    QPushButton *rButton = new QPushButton(">",vbox);
    rButton->setFixedWidth(rButton->height());
    QObject::connect(rButton,SIGNAL(pressed(void)),this,SLOT(insertChannel(void)));
    QPushButton *lButton = new QPushButton("<",vbox);
    lButton->setFixedWidth(lButton->height());
    QObject::connect(lButton,SIGNAL(pressed(void)),this,SLOT(removeChannel(void)));

    selectionBox = new QListBox(channelBox);

    sampleBox = new QVGroupBox("Sample Control",this);
    sampleBox->setInsideMargin(5);
    layout->addWidget(sampleBox);

    hbox = new QHBox(sampleBox);
    (new QLabel("Downsampling Rate:",hbox))->setFixedWidth(150);
    downsampleSpin = new QSpinBox(1,1000,1,hbox);
    QObject::connect(downsampleSpin,SIGNAL(valueChanged(int)),this,SLOT(updateDownsampleRate(int)));

    fileBox = new QVGroupBox("File Control",this);
    fileBox->setInsideMargin(5);
    layout->addWidget(fileBox);

    hbox = new QHBox(fileBox);
    hbox->setSpacing(2);
    (new QLabel("File:",hbox))->setFixedWidth(50);
    fileNameEdit = new QLineEdit(hbox);
    fileNameEdit->setReadOnly(true);
    QPushButton *fileChangeButton = new QPushButton("Change",hbox);
    QObject::connect(fileChangeButton,SIGNAL(clicked(void)),this,SLOT(changeDataFile(void)));

    hbox = new QHBox(this);
    layout->addWidget(hbox);

    startRecordButton = new QPushButton("Start Recording",hbox);
    QObject::connect(startRecordButton,SIGNAL(clicked(void)),this,SLOT(startRecordClicked(void)));
    stopRecordButton  = new QPushButton("Stop Recording",hbox);
    QObject::connect(stopRecordButton,SIGNAL(clicked(void)),this,SLOT(stopRecordClicked(void)));

    QPushButton *closeButton = new QPushButton("Close",hbox);
    QObject::connect(closeButton,SIGNAL(clicked(void)),this,SLOT(close(void)));

    resize(550,260);
    show();

    // Build initial block list
    IO::Connector::getInstance()->foreachBlock(buildBlockPtrList,&blockPtrList);
    for(std::vector<IO::Block *>::const_iterator i = blockPtrList.begin(),end = blockPtrList.end();i != end;++i)
	{
        blockList->insertItem((*i)->getName()+" "+QString::number((*i)->getID()));
	}

    // Build initial channel list
    buildChannelList();

    // Launch Recording Thread
    pthread_create(&thread,0,bounce,this);

    counter = 0;
    downsample_rate = 1;
    prev_input = 0.0;

    setActive(true);
}

DataRecorder::Panel::~Panel(void) {
    Plugin::getInstance()->removeDataRecorderPanel(this);

    setActive(false);

    DoneEvent event(fifo);
    while(RT::System::getInstance()->postEvent(&event));

    pthread_join(thread,0);

    for(RT::List<Channel>::iterator i = channels.begin(),end = channels.end();i != end;)
        delete &*(i++);
}

// ***************************************************************
// **  **  **  **  **  **  My Code  **  **  **  **  **  **  **  **
// ***************************************************************

#include <../../../fftReal/ffft/FFTReal.h>	// Several of the FFTReal files have been modified
#include <../electrode_resistance_measurement/electrode_resistance_measurement.h>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

// Parameters: Frequencies
const int SAMP_FREQ = 30.303 * 1000; 	// in Hertz
const int CUTOFF_FREQ = 6500;		// in Hz. The highest frequency considered for computing correlation coefficents
const int FREQ_LOW_CUT = 1000;		// in Hz. The lowest frequency considered.

// Parameters: RMS
const int RMS_PERIOD = 1000;		// in us
const double RMS_THRESHOLD = 0.50;	// in Volts, threshold to start recording
const int MS_IN_RMS = 10;		// in ms, time length checked to start recording
const int RMS_TERM_LEN = 200;		// in ms, time length checked to stop recording
const double RMS_TERM_THRES = .45;	// in Volts, threshold to stop recording (if under, stops)

// Parameters: FFT
const int FFT_ARRAY_SIZE = 256; 	// Must be power of 2
const int FFT_PERIOD = 1000;		// in us

// Parameters: Recording: Lengths
const int PREDATA_LENGTH = 2;		// in s; Amount of Time Before Detection to Record. Also is seconds in Circular Buffer. Make >= 1 or program will crash.
const int RECORD_H5_SECS = 5;		// in s; Min Num Seconds/ Trial Recorded
const int MAX_FILE_SECS = 25;		// in s; Max sec.s recorded, not including predata

// Parameters: Spectrogram/ Template
const char TEMP_FILE[25] = "LowerHalfTemp.Det";	// Template File (Raw Data)
const int MS_IN_SPECT = 32;		// in ms
const int SPECT_PERIOD = 1000;		// in us

// Parameters: Trigger
double corr_coeff_trig = 0.74;		// Correlation to Template needed for Detection

// Parameters: Output
const int OUTPUT_LENGTH = 32;		// in ms; output file length
const char PLAYBACK_FILE[20] = "thisSyl"; // File for Feedback
const int DEAD_TIME = 50; 		// in ms. Default: OUTPUT_LENGTH + 20

// End Parameters                      ******                      ******
// End Parameters                      ******                      ******
// End Parameters                      ******                      ******

// Circular Buffer
const int CB_ARRAY_SIZE = PREDATA_LENGTH * SAMP_FREQ;
double circularBuffer[2][CB_ARRAY_SIZE];
double predata[2][CB_ARRAY_SIZE];

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

// Physical Constants
const double PI = 3.1415926535897932;

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

	// Get File Data
	ifstream playbackFile;
	playbackFile.open(PLAYBACK_FILE);
	for (int i = 0; i != OUT_LEN_PTS; i++)
		playbackFile >> playback[i];
	playbackFile.close();

	// Get Temp Data
	ifstream tempFile;
	tempFile.open(TEMP_FILE);
	for (int i = 0; i != TD_SIZE; i++) 
		tempFile >> temp_data[i];
	tempFile.close();

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

	// Do RMS (when not recording) and Grab Pre-Data
	if (!rmsHit)
		if ((  cBcounter % (SAMP_FREQ / (1000000 / RMS_PERIOD)) == ( (SAMP_FREQ / (1000000 / RMS_PERIOD)) - 1) ))
			if ((RmsValue(POINTS_IN_RMS) >= RMS_THRESHOLD)) {
				rmsHit = true;
				for (int i = 0, j = cBcounter+1; i < CB_ARRAY_SIZE; i++, j++) {
					predata[0][i] = circularBuffer[0][j];
					predata[1][i] = circularBuffer[1][j];
					if (j == CB_ARRAY_SIZE - 1) j = -1;
				}
			}

	// Do FFT
	if (rmsHit)
		if (((cBcounter+1) % (SAMP_FREQ / (1000000 / FFT_PERIOD)) == ( (SAMP_FREQ / (1000000 / FFT_PERIOD)) - 1)))
			FftGo();

	// Compare Spectrogram with Template
	if (rmsHit)
		if (((cBcounter+2) % (SAMP_FREQ / (1000000 / SPECT_PERIOD)) == ( (SAMP_FREQ / (1000000 / SPECT_PERIOD)) - 1)))
			Compare_Spects();

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
	RMS = sqrt(RMS/ (SAMP_FREQ / (1000000 / RMS_PERIOD)) );
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

	// Zap Em? (output)
	if (correlation_coefficient >= corr_coeff_trig && deadCounter == 0) {
		isInteresting = true;
		zapEm = true;
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
        for(RT::List<Channel>::iterator i = channels.begin(),end = channels.end();i != end;++i)
            if(i->block) data[n++] = i->block->getValue(i->type,i->index);

	// Activate Command Code
	if (!writeData) Command_Code(data[0],data[1]);

	// Output Playback File
	if (zapEm) {
		if (pBcounter < OUT_LEN_PTS) 
			output(0) = playback[pBcounter];
		else if (pBcounter == OUT_LEN_PTS)
			output(0) = 0.0;
		if ((++pBcounter == OUT_LEN_PTS+1)) {
			pBcounter = 0;
			zapEm = false;
			deadCounter = 1;
		}
	}

	// Write Data
	if (writeData) {
		if (isInteresting) {
			double appendData[2];
			for (int i = 0; i != CB_ARRAY_SIZE; i++) {
				appendData[0] = predata[0][i];
				appendData[1] = predata[1][i];
				DataRecorder::Panel::AppendMe(appendData);
			}
			for (int i = 0; i != dScounter; i++) {
				appendData[0] = dataStorage[0][i];
				appendData[1] = dataStorage[1][i];
				DataRecorder::Panel::AppendMe(appendData);
			}
			cout << "  Trial " << trialNum << " Complete.\n\n"; // open RTXI from the terminal to see cout messages
			StopStart();
		}
		writeData = false;
		dScounter = 0;
	}

	// Store Data to be Written (several seconds) Later
	if (rmsHit && (dScounter < MFS_PTS)) {
		dataStorage[0][dScounter] = data[0];
		dataStorage[1][dScounter++] = data[1];
	}
    }
    // OP Code
    counter %= downsample_rate;

    // Reset, Ground Output
    if ((!recording && isRunning)) {
	Reset();
	output(0) = 0.0;
    }

    // Intialize Everything
    if (!initialized) Initialize();

    // Dead Time
    if (deadCounter != 0) deadCounter++;
    if (deadCounter >= DEAD_TIME * (SAMP_FREQ/ 1000)) deadCounter = 0;
}

void DataRecorder::Panel::receiveEvent(const Event::Object *event) {
    if(event->getName() == Event::IO_BLOCK_INSERT_EVENT) {

        IO::Block *block = reinterpret_cast<IO::Block *>(event->getParam("block"));

        blockPtrList.push_back(block);
        blockList->insertItem(block->getName()+" "+QString::number(block->getID()));
        if(blockList->count() == 1)
            buildChannelList();

    } else if(event->getName() == Event::IO_BLOCK_REMOVE_EVENT) {

        IO::Block *block = reinterpret_cast<IO::Block *>(event->getParam("block"));
        QString name = block->getName()+" "+QString::number(block->getID());

        int n = 0;
        for(;n < blockList->count() && blockList->text(n) != name;++n);
        if(n < blockList->count())
            blockList->removeItem(n);
        blockPtrList.erase(blockPtrList.begin()+n);

        for(RT::List<Channel>::iterator i = channels.begin(),end = channels.end();i != end;++i)
            if(i->block == block)
                if(recording)
                    i->block = 0;

    } else if(event->getName() == Event::START_RECORDING_EVENT) {

        StartRecordingEvent e(recording,fifo);
        RT::System::getInstance()->postEvent(&e);

    } else if(event->getName() == Event::STOP_RECORDING_EVENT) {

        StopRecordingEvent e(recording,fifo);
        RT::System::getInstance()->postEvent(&e);

    } else if(event->getName() == Event::ASYNC_DATA_EVENT) {

        AsyncDataEvent e(reinterpret_cast<double *>(event->getParam("data")),
                             *reinterpret_cast<size_t *>(event->getParam("size")),fifo);
        RT::System::getInstance()->postEvent(&e);

    }
}

void DataRecorder::Panel::receiveEventRT(const Event::Object *event) {
    if(event->getName() == Event::START_RECORDING_EVENT) {
        data_token_t token;

        recording = true;

        token.type = DataRecorder::START;
        token.size = 0;
        token.time = RT::OS::getTime();

	if (!fifo.tooBig(sizeof(token)))
	        fifo.write(&token,sizeof(token));
    } else if(event->getName() == Event::STOP_RECORDING_EVENT) {
        data_token_t token;

        recording = false;

        token.type = DataRecorder::STOP;
        token.size = 0;
        token.time = RT::OS::getTime();
	if (!fifo.tooBig(sizeof(token)))
	        fifo.write(&token,sizeof(token));
    } else if(event->getName() == Event::ASYNC_DATA_EVENT) {
        size_t size = *reinterpret_cast<size_t *>(event->getParam("size"));

        data_token_t token;

        token.type = DataRecorder::ASYNC;
        token.size = size*sizeof(double);
        token.time = RT::OS::getTime();

	if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(token.size)) {
	        fifo.write(&token,sizeof(token));
	        fifo.write(event->getParam("data"),token.size);
	}
    } else if(event->getName() == Event::WORKSPACE_PARAMETER_CHANGE_EVENT) {
        data_token_t token;

        token.type = DataRecorder::PARAM;
        token.size = sizeof(param_change_t);
        token.time = RT::OS::getTime();

        param_change_t data;
        data.id = reinterpret_cast<Settings::Object::ID>(event->getParam("object"));
        data.index = reinterpret_cast<size_t>(event->getParam("index"));
        data.step = file.idx;
        data.value = *reinterpret_cast<double *>(event->getParam("value"));

	if (!fifo.tooBig(sizeof(token)) && !fifo.tooBig(sizeof(data))) {
	        fifo.write(&token,sizeof(token));
	        fifo.write(&data,sizeof(data));
	}
    }
}

void DataRecorder::Panel::buildChannelList(void) {
    channelList->clear();

    if(!blockList->count())
        return;

    IO::Block *block = blockPtrList[blockList->currentItem()];
    IO::flags_t type;
    switch(typeList->currentItem()) {
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
          typeList->setCurrentItem(0);
          type = Workspace::INPUT;
    }

    for(size_t i = 0;i < block->getCount(type);++i)
        channelList->insertItem(block->getName(type,i));
}

void DataRecorder::Panel::changeDataFile(void) {
    QFileDialog fileDialog(this,NULL,true);
    fileDialog.setCaption("Select Data File");
    fileDialog.setMode(QFileDialog::AnyFile);

    QStringList filterList;
    filterList.push_back("HDF5 files (*.h5)");
    filterList.push_back("All files (*.*)");
    fileDialog.setFilters(filterList);

    fileDialog.exec();

    if(fileDialog.selectedFile() == "/" ||
       fileDialog.selectedFile().isNull() ||
       fileDialog.selectedFile().isEmpty())
        return;

    QString filename = fileDialog.selectedFile();

    if(!filename.lower().endsWith(QString(".h5")))
        filename += ".h5";

    OpenFileEvent event(filename,fifo);
    RT::System::getInstance()->postEvent(&event);
}

void DataRecorder::Panel::insertChannel(void) {
    if(!blockList->count() || !channelList->count())
        return;

    Channel *channel = new Channel();
    channel->block = blockPtrList[blockList->currentItem()];
    switch(typeList->currentItem()) {
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
          typeList->setCurrentItem(0);
          channel->type = Workspace::INPUT;
    }
    channel->index = channelList->currentItem();

    channel->name.sprintf("%s %ld : %s",channel->block->getName().c_str(),channel->block->getID(),
                          channel->block->getName(channel->type,channel->index).c_str());

    InsertChannelEvent event(recording,channels,channels.end(),*channel);
    if(!RT::System::getInstance()->postEvent(&event))
        selectionBox->insertItem(channel->name);
}

void DataRecorder::Panel::removeChannel(void) {
    if(!selectionBox->count())
        return;

    for(RT::List<Channel>::iterator i = channels.begin(),end = channels.end();i != end;++i)
        if(i->name == selectionBox->currentText()) {
            RemoveChannelEvent event(recording,channels,*i);
            if(!RT::System::getInstance()->postEvent(&event))
                selectionBox->removeItem(selectionBox->currentItem());
            break;
        }
}

void DataRecorder::Panel::startRecordClicked(void) {
    StartRecordingEvent event(recording,fifo);
    RT::System::getInstance()->postEvent(&event);
}

void DataRecorder::Panel::stopRecordClicked(void) {
    StopRecordingEvent event(recording,fifo);
    RT::System::getInstance()->postEvent(&event);
}

void DataRecorder::Panel::updateDownsampleRate(int r) {
    downsample_rate = r;
}

void DataRecorder::Panel::customEvent(QCustomEvent *e) {
    if(e->type() == QFileExistsEvent) {
        FileExistsEventData *data = reinterpret_cast<FileExistsEventData *>(e->data());
        data->response = QMessageBox::information(this,"File exists","The file \""+data->filename+"\" already exists.",
                                                  "Append","Overwrite","Cancel",0,2);
        data->done.wakeAll();
    } else if(e->type() == QNoFileOpenEvent) {
        QMessageBox::critical(this,"Failed to start recording",
                              "No file has been opened for writing so recording could not be started.",
                              QMessageBox::Ok,QMessageBox::NoButton);
    } else if(e->type() == QSetFileNameEditEvent) {
        SetFileNameEditEventData *data = reinterpret_cast<SetFileNameEditEventData *>(e->data());
        fileNameEdit->setText(data->filename);
        data->done.wakeAll();
    } else if(e->type() == QDisableGroupsEvent) {
        channelBox->setEnabled(false);
        sampleBox->setEnabled(false);
    } else if(e->type() == QEnableGroupsEvent) {
        channelBox->setEnabled(true);
        sampleBox->setEnabled(true);
    }
}

void DataRecorder::Panel::doDeferred(const Settings::Object::State &s) {
    for(int i = 0;i < s.loadInteger("Num Channels");++i) {
        Channel *channel;
        IO::Block *block;
        std::ostringstream str;
        str << i;

        block = dynamic_cast<IO::Block *>(Settings::Manager::getInstance()->getObject(s.loadInteger(str.str()+" ID")));
        if(!block) continue;

        channel = new Channel();

        channel->block = block;
        channel->type = s.loadInteger(str.str()+" type");
        channel->index = s.loadInteger(str.str()+" index");
        channel->name.sprintf("%s %ld : %s",channel->block->getName().c_str(),channel->block->getID(),
                              channel->block->getName(channel->type,channel->index).c_str());

        channels.insert(channels.end(),*channel);
        selectionBox->insertItem(channel->name);
    }
}

void DataRecorder::Panel::doLoad(const Settings::Object::State &s) {
    if(s.loadInteger("Maximized"))
        showMaximized();
    else if(s.loadInteger("Minimized"))
        showMinimized();

    downsampleSpin->setValue(s.loadInteger("Downsample"));

    resize(s.loadInteger("W"),
           s.loadInteger("H"));
    parentWidget()->move(s.loadInteger("X"),
                         s.loadInteger("Y"));
}

void DataRecorder::Panel::doSave(Settings::Object::State &s) const {
    if(isMaximized())
        s.saveInteger("Maximized",1);
    else if(isMinimized())
        s.saveInteger("Minimized",1);

    QPoint pos = parentWidget()->pos();
    s.saveInteger("X",pos.x());
    s.saveInteger("Y",pos.y());
    s.saveInteger("W",width());
    s.saveInteger("H",height());

    s.saveInteger("Downsample",downsampleSpin->value());

    s.saveInteger("Num Channels",channels.size());
    size_t n = 0;
    for(RT::List<Channel>::const_iterator i = channels.begin(),end = channels.end();i != end;++i) {
        std::ostringstream str;
        str << n++;

        s.saveInteger(str.str()+" ID",i->block->getID());
        s.saveInteger(str.str()+" type",i->type);
        s.saveInteger(str.str()+" index",i->index);
    }
}

void *DataRecorder::Panel::bounce(void *param) {
    Panel *that = reinterpret_cast<Panel *>(param);
    if(that) that->processData();
    return 0;
}

void DataRecorder::Panel::processData(void) {
    enum {
        CLOSED,
        OPENED,
        RECORD,
    } state = CLOSED;

    data_token_t token;

    for(;;) {

        fifo.read(&token,sizeof(token));

        if(token.type == SYNC) {

            if(state == RECORD) {
                double data[token.size/sizeof(double)];
                fifo.read(data,token.size);
                H5PTappend(file.cdata,1,data);

                ++file.idx;
            }

        } else if(token.type == ASYNC) {

            if(state == RECORD) {

                double data[token.size/sizeof(double)];
                fifo.read(data,token.size);

                if(data) {
                    hsize_t array_size[] = { token.size/sizeof(double) };
                    hid_t array_space = H5Screate_simple(1,array_size,array_size);
                    hid_t array_type = H5Tarray_create(H5T_IEEE_F64LE,1,array_size);

                    QString data_name = QString::number(static_cast<unsigned long long>(token.time));

                    hid_t adata = H5Dcreate(file.adata,data_name.latin1(),array_type,array_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
                    H5Dwrite(adata,array_type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);

                    H5Dclose(adata);
                    H5Tclose(array_type);
                    H5Sclose(array_space);
                }
            }

        } else if(token.type == OPEN) {

            if(state == RECORD)
                stopRecording(token.time);

            if(state != CLOSED)
                closeFile();

            char filename_string[token.size];
            fifo.read(filename_string,token.size);
            QString filename = filename_string;

            if(openFile(filename))
                state = CLOSED;
            else
                state = OPENED;

        } else if(token.type == CLOSE) {

            if(state == RECORD)
                stopRecording(RT::OS::getTime());

            if(state != CLOSED)
                closeFile();

            state = CLOSED;

        } else if(token.type == START) {

            if(state == CLOSED) {
                QCustomEvent *event = new QCustomEvent(QNoFileOpenEvent);
                QApplication::postEvent(this,event);
            } else if(state == OPENED) {
                startRecording(token.time);
                state = RECORD;
            }

        } else if(token.type == STOP) {

            if(state == RECORD) {
                stopRecording(token.time);
                state = OPENED;
            }

        } else if(token.type == DONE) {

            if(state == RECORD)
                stopRecording(token.time,true);

            if(state != CLOSED)
                closeFile(true);

            break;
        } else if(token.type == PARAM) {
            param_change_t data;
            fifo.read(&data,sizeof(data));

            IO::Block *block = dynamic_cast<IO::Block *>(Settings::Manager::getInstance()->getObject(data.id));

            if(block && state == RECORD) {
                param_hdf_t param = {
                    data.step,
                    data.value,
                };

                hid_t param_type;
                param_type = H5Tcreate(H5T_COMPOUND,sizeof(param_hdf_t));
                H5Tinsert(param_type,"index",HOFFSET(param_hdf_t,index),H5T_STD_I64LE);
                H5Tinsert(param_type,"value",HOFFSET(param_hdf_t,value),H5T_IEEE_F64LE);

                QString parameter_name = QString::number(block->getID())+" "+block->getName()+" : "+block->getName(Workspace::PARAMETER,data.index);

                hid_t data = H5PTopen(file.pdata,parameter_name.latin1());
                H5PTappend(data,1,&param);
                H5PTclose(data);

                H5Tclose(param_type);
            }
        }

    }
}

int DataRecorder::Panel::openFile(QString &filename) {

#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread)) {
        ERROR_MSG("DataRecorder::Panel::openFile : called by invalid thread\n");
        PRINT_BACKTRACE();
    }
#endif

    if(QFile::exists(filename)) {
        QCustomEvent *event = new QCustomEvent(QFileExistsEvent);
        FileExistsEventData data;

        event->setData(&data);
        data.filename = filename;

        QApplication::postEvent(this,event);
        data.done.wait();

        if(data.response == 0)
            file.id = H5Fopen(filename.latin1(),H5F_ACC_RDWR,H5P_DEFAULT);
        else if(data.response == 1)
            file.id = H5Fcreate(filename.latin1(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        else
            return -1;
    } else
        file.id = H5Fcreate(filename.latin1(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);

    if(file.id < 0) {
        H5E_type_t error_type;
        size_t error_size;
        error_size = H5Eget_msg(file.id,&error_type,NULL,0);

        {
            char error_msg[error_size+1];
            H5Eget_msg(file.id,&error_type,error_msg,error_size);
            error_msg[error_size] = 0;
            H5Eclear(file.id);

            ERROR_MSG("DataRecorder::Panel::processData : failed to open \"%s\" for writing with error : %s\n",filename.latin1(),error_msg);
            return -1;
        }
    }

    QCustomEvent *event = new QCustomEvent(QSetFileNameEditEvent);
    SetFileNameEditEventData data;

    event->setData(&data);
    data.filename = filename;

    QApplication::postEvent(this,event);
    data.done.wait();

    return 0;
}

void DataRecorder::Panel::closeFile(bool shutdown) {

#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread)) {
        ERROR_MSG("DataRecorder::Panel::closeFile : called by invalid thread\n");
        PRINT_BACKTRACE();
    }
#endif

    H5Fclose(file.id);

    if(!shutdown) {
        QCustomEvent *event = new QCustomEvent(QSetFileNameEditEvent);
        SetFileNameEditEventData data;

        event->setData(&data);
        data.filename = "";

        QApplication::postEvent(this,event);
        data.done.wait();
    }
}

int DataRecorder::Panel::startRecording(long long timestamp) {

#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread)) {
        ERROR_MSG("DataRecorder::Panel::startRecording : called by invalid thread\n");
        PRINT_BACKTRACE();
    }
#endif

    size_t trial_num;
    QString trial_name;

    H5Eset_auto(H5E_DEFAULT,NULL,NULL);


    for(trial_num=1;;++trial_num) {
	trial_name = "/Trial"+QString::number(trial_num);
	file.trial = H5Gopen(file.id,trial_name.latin1(),H5P_DEFAULT);

	if((file.trial < 0)) {
	    H5Eclear(H5E_DEFAULT);
	    break;
	} else
	    H5Gclose(file.trial);
    }
	trialNum = trial_num;


    file.trial = H5Gcreate(file.id,trial_name.latin1(),H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    file.pdata = H5Gcreate(file.trial,"Parameters",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    file.adata = H5Gcreate(file.trial,"Asyncronous Data",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    file.sdata = H5Gcreate(file.trial,"Syncronous Data",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t string_type = H5Tcopy(H5T_C_S1);
    size_t string_size = 512;
    H5Tset_size(string_type,string_size);
    hid_t data;

    long long period = RT::System::getInstance()->getPeriod();
    data = H5Dcreate(file.trial,"Period (ns)",H5T_STD_U64LE,scalar_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(data,H5T_STD_U64LE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&period);
    H5Dclose(data);

    long long downsample = downsample_rate;
    data = H5Dcreate(file.trial,"Downsampling Rate",H5T_STD_U64LE,scalar_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(data,H5T_STD_U64LE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&downsample);
    H5Dclose(data);

    data = H5Dcreate(file.trial,"Timestamp Start (ns)",H5T_STD_U64LE,scalar_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(data,H5T_STD_U64LE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&timestamp);
    H5Dclose(data);

    data = H5Dcreate(file.trial,"Date",string_type,scalar_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(data,string_type,H5S_ALL,H5S_ALL,H5P_DEFAULT,QDateTime::currentDateTime().toString(Qt::ISODate).latin1());
    H5Dclose(data);

    hid_t param_type;
    param_type = H5Tcreate(H5T_COMPOUND,sizeof(param_hdf_t));
    H5Tinsert(param_type,"index",HOFFSET(param_hdf_t,index),H5T_STD_I64LE);
    H5Tinsert(param_type,"value",HOFFSET(param_hdf_t,value),H5T_IEEE_F64LE);

    for(RT::List<Channel>::iterator i = channels.begin(), end = channels.end();i != end;++i) {
        IO::Block *block = i->block;
        for(size_t j = 0;j < block->getCount(Workspace::PARAMETER);++j) {
            QString parameter_name = QString::number(block->getID())+" "+block->getName()+" : "+block->getName(Workspace::PARAMETER,j);
            data = H5PTcreate_fl(file.pdata,parameter_name.latin1(),param_type,sizeof(param_hdf_t),-1);
            struct param_hdf_t value = {
                0,
                block->getValue(Workspace::PARAMETER,j),
            };
            H5PTappend(data,1,&value);
            H5PTclose(data);
        }
        for(size_t j = 0;j < block->getCount(Workspace::COMMENT);++j) {
            QString comment_name = QString::number(block->getID())+" "+block->getName()+" : "+block->getName(Workspace::COMMENT,j);
            hsize_t dims = dynamic_cast<Workspace::Instance *>(block)->getValueString(Workspace::COMMENT,j).size()+1;
            hid_t comment_space = H5Screate_simple(1,&dims,&dims);
            data = H5Dcreate(file.pdata,comment_name.latin1(),H5T_C_S1,comment_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
            H5Dwrite(data,H5T_C_S1,H5S_ALL,H5S_ALL,H5P_DEFAULT,dynamic_cast<Workspace::Instance *>(block)->getValueString(Workspace::COMMENT,j).c_str());
            H5Dclose(data);
        }
    }

    H5Tclose(param_type);

    size_t count = 0;
    for(RT::List<Channel>::iterator i = channels.begin(), end = channels.end();i != end;++i) {
        QString channel_name = "Channel " + QString::number(++count) + " Name";
        hid_t data = H5Dcreate(file.sdata,channel_name.latin1(),string_type,scalar_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(data,string_type,H5S_ALL,H5S_ALL,H5P_DEFAULT,i->name.latin1());
        H5Dclose(data);
    }

    H5Tclose(string_type);
    H5Sclose(scalar_space);

    if(channels.size()) {
        hsize_t array_size[] = { channels.size() };
        hid_t array_type = H5Tarray_create(H5T_IEEE_F64LE,1,array_size);
        file.cdata = H5PTcreate_fl(file.sdata,"Channel Data",array_type,(hsize_t)64,1);
        H5Tclose(array_type);
    }

    file.idx = 0;

    QCustomEvent *event = new QCustomEvent(QDisableGroupsEvent);
    QApplication::postEvent(this,event);

    return 0;
}

void DataRecorder::Panel::stopRecording(long long timestamp,bool shutdown) {

#ifdef DEBUG
    if(!pthread_equal(pthread_self(),thread)) {
        ERROR_MSG("DataRecorder::Panel::stopRecording : called by invalid thread\n");
        PRINT_BACKTRACE();
    }
#endif

    hid_t scalar_space = H5Screate(H5S_SCALAR);
    hid_t data = H5Dcreate(file.trial,"Timestamp Stop (ns)",H5T_STD_U64LE,scalar_space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(data,H5T_STD_U64LE,H5S_ALL,H5S_ALL,H5P_DEFAULT,&timestamp);
    H5Dclose(data);
    H5Sclose(scalar_space);

    H5PTclose(file.cdata);
    H5Gclose(file.sdata);
    H5Gclose(file.pdata);
    H5Gclose(file.adata);
    H5Gclose(file.trial);

    H5Fflush(file.id,H5F_SCOPE_LOCAL);
    void *file_handle;
    H5Fget_vfd_handle(file.id,H5P_DEFAULT,&file_handle);
    if(fsync(*static_cast<int *>(file_handle))) {
        DEBUG_MSG("DataRecorder::Panel::stopRecording : fsync failed, running sync\n");
        sync();
    }

    if(!shutdown) {
        QCustomEvent *event = new QCustomEvent(QEnableGroupsEvent);
        QApplication::postEvent(this,event);
    }
}

extern "C" Plugin::Object *createRTXIPlugin(void *) {
    return DataRecorder::Plugin::getInstance();
}

DataRecorder::Plugin::Plugin(void) {
    menuID = MainWindow::getInstance()->createControlMenuItem("Data Recorder",this,SLOT(createDataRecorderPanel(void)));
}

DataRecorder::Plugin::~Plugin(void) {
    MainWindow::getInstance()->removeControlMenuItem(menuID);
    while(panelList.size())
        delete panelList.front();
    instance = 0;
}

DataRecorder::Panel *DataRecorder::Plugin::createDataRecorderPanel(void) {
    Panel *panel = new Panel(MainWindow::getInstance()->centralWidget());
    panelList.push_back(panel);
    return panel;
}

void DataRecorder::Plugin::removeDataRecorderPanel(DataRecorder::Panel *panel) {
    panelList.remove(panel);
}

void DataRecorder::Plugin::doDeferred(const Settings::Object::State &s) {
    size_t i = 0;
    for(std::list<Panel *>::iterator j = panelList.begin(),end = panelList.end();j != end;++j)
        (*j)->deferred(s.loadState(QString::number(i++)));
}

void DataRecorder::Plugin::doLoad(const Settings::Object::State &s) {
    for(size_t i=0;i < static_cast<size_t>(s.loadInteger("Num Panels"));++i) {
        Panel *panel = new Panel(MainWindow::getInstance()->centralWidget());
        panelList.push_back(panel);
        panel->load(s.loadState(QString::number(i)));
    }
}

void DataRecorder::Plugin::doSave(Settings::Object::State &s) const {
    s.saveInteger("Num Panels",panelList.size());
    size_t n = 0;
    for(std::list<Panel *>::const_iterator i = panelList.begin(),end = panelList.end();i != end;++i)
        s.saveState(QString::number(n++),(*i)->save());
}

static Mutex mutex;
DataRecorder::Plugin *DataRecorder::Plugin::instance = 0;

DataRecorder::Plugin *DataRecorder::Plugin::getInstance(void) {
    if(instance)
        return instance;

    /*************************************************************************
     * Seems like alot of hoops to jump through, but allocation isn't        *
     *   thread-safe. So effort must be taken to ensure mutual exclusion.    *
     *************************************************************************/

    Mutex::Locker lock(&::mutex);
    if(!instance)
        instance = new Plugin();

    return instance;
}
