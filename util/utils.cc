// $Id: utils.cc,v 1.1 2011/02/17 13:26:13 mschrode Exp $

#include "ConfigParser.h"
#include "FileOps.h"
#include "HistOps.h"
#include "LabelFactory.h"
#include "StyleSettings.h"
#include "utils.h"


std::vector<TH1*> util::LabelFactory::hDummies_;
unsigned int util::HistOps::COUNT_ = 0;





