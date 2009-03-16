/*
 * Copyright (C) 2006 Weill Medical College of Cornell University
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

#include <connector.h>
#include <main_window.h>
#include <mutex.h>
#include <qcombobox.h>
#include <qgroupbox.h>
#include <qhbox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlistbox.h>
#include <qpushbutton.h>

struct block_list_info_t {
    QComboBox *blockList0;
    QComboBox *blockList1;
    std::vector<IO::Block *> *blocks;
};

static void buildBlockList(IO::Block *block,void *arg) {
    block_list_info_t *info = static_cast<block_list_info_t *>(arg);
    info->blockList0->insertItem(block->getName()+QString(" ")+QString::number(block->getID()));
    info->blockList1->insertItem(block->getName()+QString(" ")+QString::number(block->getID()));
    info->blocks->push_back(block);
}

Connector::Panel::Panel(QWidget *parent)
    : QWidget(parent) {

    setCaption("Connector Panel");

    QGridLayout *layout = new QGridLayout(this,2,3);

    QGroupBox *gbox0 = new QGroupBox("Output Block",this);
    layout->addWidget(gbox0,0,0);
    QBoxLayout *boxLayout0 = new QVBoxLayout(gbox0,5,5);
    boxLayout0->setAutoAdd(true);
    boxLayout0->addSpacing(18);

    QHBox *hbox0 = new QHBox(gbox0);
    (new QLabel("Block: ",hbox0))->setFixedWidth(75);
    outputBlock = new QComboBox(hbox0);
    QObject::connect(outputBlock,SIGNAL(activated(int)),this,SLOT(buildOutputChannelList(void)));

    QHBox *hbox1 = new QHBox(gbox0);
    (new QLabel("Channel: ",hbox1))->setFixedWidth(75);
    outputChannel = new QComboBox(hbox1);
    QObject::connect(outputChannel,SIGNAL(activated(int)),this,SLOT(updateConnectionButton(void)));

    connectionButton = new QPushButton("==>",this);
    layout->addWidget(connectionButton,0,1);
    connectionButton->setToggleButton(true);
    QObject::connect(connectionButton,SIGNAL(toggled(bool)),this,SLOT(toggleConnection(bool)));

    QGroupBox *gbox1 = new QGroupBox("Input Block",this);
    layout->addWidget(gbox1,0,2);
    QBoxLayout *boxLayout1 = new QVBoxLayout(gbox1,5,5);
    boxLayout1->setAutoAdd(true);
    boxLayout1->addSpacing(18);

    QHBox *hbox2 = new QHBox(gbox1);
    (new QLabel("Block: ",hbox2))->setFixedWidth(75);
    inputBlock = new QComboBox(hbox2);
    QObject::connect(inputBlock,SIGNAL(activated(int)),this,SLOT(buildInputChannelList(void)));

    QHBox *hbox3 = new QHBox(gbox1);
    (new QLabel("Channel: ",hbox3))->setFixedWidth(75);
    inputChannel = new QComboBox(hbox3);
    QObject::connect(inputChannel,SIGNAL(activated(int)),this,SLOT(updateConnectionButton(void)));

    QGroupBox *gbox2 = new QGroupBox("Connections",this);
    layout->addMultiCellWidget(gbox2,1,1,0,2);
    QVBoxLayout *boxLayout2 = new QVBoxLayout(gbox2,5,5);
    boxLayout2->setAutoAdd(true);
    boxLayout2->addSpacing(18);

    connectionBox = new QListBox(gbox2);
    QObject::connect(connectionBox,SIGNAL(selected(const QString &)),this,SLOT(highlightConnectionBox(const QString &)));

    block_list_info_t info = { inputBlock, outputBlock, &blocks };
    IO::Connector::foreachBlock(::buildBlockList,&info);

    if(blocks.size() >= 1) {
        buildInputChannelList();
        buildOutputChannelList();
    }

    IO::Connector::foreachConnection(&buildConnectionList,&links);
    for(size_t i = 0, iend = links.size();i < iend;++i) {
        connectionBox->insertItem(QString::number(links[i].src->getID())+" "+links[i].src->getName()+" : "+QString::number(links[i].src_idx)+" "+links[i].src->getName(IO::OUTPUT,links[i].src_idx)+" ==> "+
                              QString::number(links[i].dest->getID())+" "+links[i].dest->getName()+" : "+QString::number(links[i].dest_idx)+" "+links[i].dest->getName(IO::INPUT,links[i].dest_idx));
    }
}

Connector::Panel::~Panel(void) {}

void Connector::Panel::receiveEvent(const Event::Object *event) {
    if(event->getName() == Event::IO_BLOCK_INSERT_EVENT) {
        IO::Block *block = reinterpret_cast<IO::Block *>(event->getParam("block"));

        inputBlock->insertItem(block->getName()+QString(" ")+QString::number(block->getID()));
        outputBlock->insertItem(block->getName()+QString(" ")+QString::number(block->getID()));
        blocks.push_back(block);

        if(blocks.size() == 1) {
            buildInputChannelList();
            buildOutputChannelList();
        }
    } else if(event->getName() == Event::IO_BLOCK_REMOVE_EVENT) {
        IO::Block *block = reinterpret_cast<IO::Block *>(event->getParam("block"));

        size_t index;
        for(index = 0;index < blocks.size() && blocks[index] != block;++index);
        if(index >= blocks.size())
            return;

        size_t current0 = inputBlock->currentItem();
        size_t current1 = outputBlock->currentItem();

        inputBlock->removeItem(index);
        outputBlock->removeItem(index);
        blocks.erase(blocks.begin()+index);

        if(current0 == index) {
            inputBlock->setCurrentItem(0);
            buildInputChannelList();
        }
        if(current1 == index) {
            outputBlock->setCurrentItem(0);
            buildOutputChannelList();
        }
    } else if(event->getName() == Event::IO_LINK_INSERT_EVENT) {
        IO::Block *src = reinterpret_cast<IO::Block *>(event->getParam("src"));
        size_t src_idx = *reinterpret_cast<size_t *>(event->getParam("src_num"));
        IO::Block *dest = reinterpret_cast<IO::Block *>(event->getParam("dest"));
        size_t dest_idx = *reinterpret_cast<size_t *>(event->getParam("dest_num"));

        connectionBox->insertItem(QString::number(src->getID())+" "+src->getName()+" : "+
                                  QString::number(src_idx)+" "+src->getName(IO::OUTPUT,src_idx)+" ==> "+
                                  QString::number(dest->getID())+" "+dest->getName()+" : "+
                                  QString::number(dest_idx)+" "+dest->getName(IO::INPUT,dest_idx));
    } else if(event->getName() == Event::IO_LINK_REMOVE_EVENT) {
        IO::Block *src = reinterpret_cast<IO::Block *>(event->getParam("src"));
        size_t src_idx = *reinterpret_cast<size_t *>(event->getParam("src_num"));
        IO::Block *dest = reinterpret_cast<IO::Block *>(event->getParam("dest"));
        size_t dest_idx = *reinterpret_cast<size_t *>(event->getParam("dest_num"));

        QString link_name = QString::number(src->getID())+" "+src->getName()+" : "+
                            QString::number(src_idx)+" "+src->getName(IO::OUTPUT,src_idx)+" ==> "+
                            QString::number(dest->getID())+" "+dest->getName()+" : "+
                            QString::number(dest_idx)+" "+dest->getName(IO::INPUT,dest_idx);

        size_t index;
        for(index=0;index < connectionBox->count() && connectionBox->text(index) != link_name;++index);
        if(index >= connectionBox->count())
            ERROR_MSG("Connector::Panel::receiveEvent : removing non-existant link.\n");
        else
            connectionBox->removeItem(index);
    }
}

void Connector::Panel::buildInputChannelList(void) {
    inputChannel->clear();
    if(!inputBlock->count())
        return;

    IO::Block *block = blocks[inputBlock->currentItem()];

    for(size_t i = 0;i < block->getCount(IO::INPUT);++i)
        inputChannel->insertItem(block->getName(IO::INPUT,i));

    updateConnectionButton();
}

void Connector::Panel::buildOutputChannelList(void) {
    outputChannel->clear();
    if(!outputBlock->count())
        return;

    IO::Block *block = blocks[outputBlock->currentItem()];

    for(size_t i = 0;i < block->getCount(IO::OUTPUT);++i)
        outputChannel->insertItem(block->getName(IO::OUTPUT,i));

    updateConnectionButton();
}

void Connector::Panel::highlightConnectionBox(const QString &s) {
    QString selection = s;
    QString substr;
    int sep;

    sep = selection.find(' ');
    substr = selection.left(sep);
    selection = selection.right(selection.length()-sep-1);
    Settings::Object::ID src_id = substr.toULong();

    sep = selection.find(':');
    selection = selection.right(selection.length()-sep-2);

    sep = selection.find(' ');
    substr = selection.left(sep);
    selection = selection.right(selection.length()-sep-1);
    size_t src_idx = substr.toULong();

    sep = selection.find("==>");
    selection = selection.right(selection.length()-sep-4);

    sep = selection.find(' ');
    substr = selection.left(sep);
    selection = selection.right(selection.length()-sep-1);
    Settings::Object::ID dest_id = substr.toULong();

    sep = selection.find(':');
    selection = selection.right(selection.length()-sep-2);

    sep = selection.find(' ');
    substr = selection.left(sep);
    selection = selection.right(selection.length()-sep-1);
    size_t dest_idx = substr.toULong();

    IO::Block *src = dynamic_cast<IO::Block *>(Settings::Manager::getObject(src_id));

    size_t index;
    for(index = 0;index < blocks.size() && blocks[index] != src;++index);
    if(index >= blocks.size())
        ERROR_MSG("Connector::Panel::highlightConnectionBox : highlighted source does not exist.\n");

    outputBlock->setCurrentItem(index);
    buildOutputChannelList();
    outputChannel->setCurrentItem(src_idx);

    IO::Block *dest = dynamic_cast<IO::Block *>(Settings::Manager::getObject(dest_id));
    for(index = 0;index < blocks.size() && blocks[index] != dest;++index);
    if(index >= blocks.size())
        ERROR_MSG("Connector::Panel::highlightConnectionBox : highlighted destination does not exist.\n");

    inputBlock->setCurrentItem(index);
    buildInputChannelList();
    inputChannel->setCurrentItem(dest_idx);
}

void Connector::Panel::toggleConnection(bool on) {
    IO::Block *src = blocks[outputBlock->currentItem()];
    IO::Block *dest = blocks[inputBlock->currentItem()];
    size_t src_num = outputChannel->currentItem();
    size_t dest_num = inputChannel->currentItem();

    if(IO::Connector::connected(src,src_num,dest,dest_num) == on)
        return;

    if(on)
        IO::Connector::connect(src,src_num,dest,dest_num);
    else
        IO::Connector::disconnect(src,src_num,dest,dest_num);
}

void Connector::Panel::updateConnectionButton(void) {
    IO::Block *src = blocks[outputBlock->currentItem()];
    IO::Block *dest = blocks[inputBlock->currentItem()];
    size_t src_num = outputChannel->currentItem();
    size_t dest_num = inputChannel->currentItem();

    connectionButton->setOn(IO::Connector::connected(src,src_num,dest,dest_num));
}

void Connector::Panel::buildConnectionList(IO::Block *src,size_t src_num,IO::Block *dest,size_t dest_num,void *arg) {
    std::vector<link_t> &list = *reinterpret_cast<std::vector<link_t> *>(arg);

    link_t link = {
        src,
        src_num,
        dest,
        dest_num,
    };

    list.push_back(link);
}

extern "C" Plugin::Object *createRTXIPlugin(void *) {
    return Connector::Plugin::getInstance();
}

Connector::Plugin::Plugin(void)
    : panel(0) {
    menuID = MainWindow::createControlMenuItem("Connector",this,SLOT(showConnectorPanel(void)));
}

Connector::Plugin::~Plugin(void) {
    MainWindow::removeControlMenuItem(menuID);
    if(panel)
        delete panel;
    instance = 0;
}

void Connector::Plugin::showConnectorPanel(void) {
    if(!panel)
        panel = new Panel(MainWindow::getInstance()->centralWidget());
    panel->show();
}

static Mutex mutex;
Connector::Plugin *Connector::Plugin::instance = 0;

Connector::Plugin *Connector::Plugin::getInstance(void) {
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
