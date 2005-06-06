#include "SiconosException.h"
using namespace std;

SiconosException::SiconosException(): reportMsg("Siconos Exception")
{}

SiconosException::SiconosException(const string& report): reportMsg(report)
{}

SiconosException::~SiconosException() {}
