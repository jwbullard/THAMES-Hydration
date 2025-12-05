/**
@file Exceptions.cc
@brief Definition of methods for the various exception classes.
*/

#include "Exceptions.h"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

// Definition the EOBException class
EOBException::EOBException() {
  classname_ = "";
  functionname_ = "";
  arrayname_ = "";
  sizelimit_ = 0;
  indx_ = 0;
}

EOBException::EOBException(const string &cname, const string &fileName,
                           const string &arname, const int sl,
                           const unsigned int id) {
  classname_ = cname;
  functionname_ = fileName;
  arrayname_ = arname;
  sizelimit_ = sl;
  indx_ = id;
}

void EOBException::printException() {
  // std::clog << endl << "EOB Exception Thrown:" << endl;
  // std::clog << "    Details: " << endl;
  // std::clog << "        Offending Function " << classname_ << "::" <<
  // functionname_
  //     << endl;
  cerr << endl << "EOB Exception Thrown:" << endl;
  cerr << "    Details: " << endl;
  cerr << "        Offending Function " << classname_ << "::" << functionname_
       << endl;
  if (indx_ == 0) {
    // std::clog << "        Array: " << arrayname_ << endl;
    cerr << "        Array: " << arrayname_ << endl;
  } else {
    // std::clog << "        Array: " << arrayname_ << " contains " <<
    // sizelimit_; std::clog << " elements, but tried to access element " <<
    // indx_ << endl;
    cerr << "        Array: " << arrayname_ << " contains " << sizelimit_;
    cerr << " elements, but tried to access element " << indx_ << endl;
  }
  return;
}
// End of the EOBException class

// Definition the FileException class
FileException::FileException() {
  classname_ = "";
  functionname_ = "";
  filename_ = "";
  extype_ = "";
}

FileException::FileException(const string &cname, const string &fileName,
                             const string &filename, const string &extype) {
  classname_ = cname;
  functionname_ = fileName;
  filename_ = filename;
  extype_ = extype;
}

void FileException::printException() {
  // std::clog << endl << "File Exception Thrown:" << endl;
  // std::clog << "    Details: " << endl;
  // std::clog << "        Offending Function " << classname_ << "::" <<
  // functionname_
  //      << endl;
  // std::clog << "        File: " << filename_ << ", Problem:" << extype_ <<
  // endl;
  cerr << endl << "File Exception Thrown:" << endl;
  cerr << "    Details: " << endl;
  cerr << "        Offending Function " << classname_ << "::" << functionname_
       << endl;
  cerr << "        File: " << filename_ << ", Problem: " << extype_ << endl;
  return;
}
// End of the FileException class

// Definition the FloatException class
FloatException::FloatException() {
  classname_ = "";
  functionname_ = "";
  description_ = "";
}

FloatException::FloatException(const string &cname, const string &fileName,
                               const string &strd) {
  classname_ = cname;
  functionname_ = fileName;
  description_ = strd;
}

void FloatException::printException() {
  // std::clog << endl << "Floating Point Exception Thrown:" << endl;
  // std::clog << "    Details: " << endl;
  // std::clog << "        Offending Function " << classname_ << "::" <<
  // functionname_
  //      << endl;
  // std::clog << "        Description: " << description_ << endl;
  cerr << endl << "Floating Point Exception Thrown:" << endl;
  cerr << "    Details: " << endl;
  cerr << "        Offending Function " << classname_ << "::" << functionname_
       << endl;
  cerr << "        Description: " << description_ << endl;
  return;
}
// End of the FloatException class

// Definition the HandleException class
HandleException::HandleException() {
  classname_ = "";
  functionname_ = "";
  handle_ = "";
  description_ = "";
}

HandleException::HandleException(const string &cname, const string &fileName,
                                 const string &handle, const string &strd) {
  classname_ = cname;
  functionname_ = fileName;
  handle_ = handle;
  description_ = strd;
}

void HandleException::printException() {
  // std::clog << endl << "Handle Exception Thrown:" << endl;
  // std::clog << "    Details: " << endl;
  // std::clog << "        Offending Function " << classname_ << "::" <<
  // functionname_
  //      << endl;
  // std::clog << "        Description: " << description_ << endl;
  // std::clog << "             Handle: " << handle_ << endl;
  cerr << endl << "Floating Point Exception Thrown:" << endl;
  cerr << "    Details: " << endl;
  cerr << "        Offending Function " << classname_ << "::" << functionname_
       << endl;
  cerr << "        Description: " << description_ << endl;
  cerr << "             Handle: " << handle_ << endl;
  return;
}
// End of the HandleException class

// Definition the GEMException class
GEMException::GEMException() {
  classname_ = "";
  functionname_ = "";
  description_ = "";
}

/**
@brief Overloaded constructor that is typically invoked by THAMES.

@param cname is the class name where the exception was thrown
@param fileName is the method name where the exception was thrown
@param strd is the description of the exception
*/
GEMException::GEMException(const string &cname, const string &fileName,
                           const string &strd) {
  classname_ = cname;
  functionname_ = fileName;
  description_ = strd;
}

void GEMException::printException() {
  cerr << endl << "GEM Exception Thrown:" << endl;
  cerr << "    Details: " << endl;
  cerr << "        Offending Function " << classname_ << "::" << functionname_
       << endl;
  cerr << "        " << description_ << endl;
  return;
}
// End of GEMException class

// Definition the MicrostructureException class
MicrostructureException::MicrostructureException() {
  classname_ = "";
  functionname_ = "";
  description_ = "";
  excp_ = true;
}

/**
  @brief Overloaded constructor that is typically invoked by THAMES.

  @param cname is the class name where the exception was thrown
  @param fileName is the method name where the exception was thrown
  @param strd is the description of the exception
*/
MicrostructureException::MicrostructureException(const string &cname,
                                                 const string &fileName,
                                                 const string &strd) {
  classname_ = cname;
  functionname_ = fileName;
  description_ = strd;
}

MicrostructureException::MicrostructureException(const string &cname,
                                                 const string &fileName,
                                                 const string &strd,
                                                 bool excp) {
  classname_ = cname;
  functionname_ = fileName;
  description_ = strd;
  excp_ = excp;
}

void MicrostructureException::printException() {
  // bool excp1_ = true;
  if (excp_) {
    cerr << endl << "Microstructure Exception Thrown:" << endl;
    cerr << "    Details: " << endl;
    cerr << "        Offending Function " << classname_ << "::" << functionname_
         << endl;
    cerr << "        Problem: " << description_ << endl;
  } else {
    // std::clog << endl << "Microstructure Exception Thrown:" << endl;
    // std::clog << "    Details: " << endl;
    // std::clog << "        From Function " << classname_ << "::" <<
    // functionname_
    //      << endl;
    // std::clog << "        reason: " << description_ << endl;
    cerr << endl << "Microstructure Exception Thrown:" << endl;
    cerr << "    Details: " << endl;
    cerr << "        From Function " << classname_ << "::" << functionname_
         << endl;
    cerr << "        reason: " << description_ << endl;
  }

  return;
}
// End of MicrostructureException class

// Definition the DataException class
DataException::DataException() {
  classname_ = "";
  functionname_ = "";
  description_ = "";
}

DataException::DataException(const string &cname, const string &functionName,
                             const string &strd) {
  classname_ = cname;
  functionname_ = functionName;
  description_ = strd;
}

void DataException::printException() {
  // std::clog << endl << "Data Exception Thrown:" << endl;
  // std::clog << "    Details: " << endl;
  // std::clog << "        Offending Function " << classname_ << "::" <<
  // functionname_
  //      << endl;
  // std::clog << "        Problem:" << description_ << endl;
  cerr << endl << "Data Exception Thrown:" << endl;
  cerr << "    Details: " << endl;
  cerr << "        Offending Function " << classname_ << "::" << functionname_
       << endl;
  cerr << "        Problem: " << description_ << endl;
  return;
}
// End of DataException class
