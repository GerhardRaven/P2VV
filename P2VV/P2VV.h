/******************************************************************************
** P2VV:                                                                     **
** common P2VV tools                                                         **
**                                                                           **
** authors:                                                                  **
**   JvL, Jeroen van Leerdam, Nikhef, j.van.leerdam@nikhef.nl                **
**   WH,  Wouter Hulsbergen,  Nikhef                                         **
**                                                                           **
******************************************************************************/

#ifndef P2VV_H
#define P2VV_H
#include <iostream>
#include <string>

class ProgressDisplay {
public:
  // ProgressDisplay: os is hint;
  // implementation may ignore, particularly in embedded systems
  explicit ProgressDisplay(unsigned long expCount,
                           std::ostream & os = std::cout,
                           const std::string& s1 = "\n", // leading strings
                           const std::string& s2 = "",
                           const std::string& s3 = "") :
    _os(os), _s1(s1), _s2(s2), _s3(s3) {restart(expCount);}

  // restart effects: display appropriate scale
  // postconditions: count()==0, expCount()==expCount
  void restart(unsigned long expCount);

  //  increment effects: display appropriate progress tic if needed
  //  postconditions: count()== original count() + increment
  //  returns: count().
  unsigned long  increment(unsigned long incr);

  unsigned long  operator+=(unsigned long incr) {return increment(incr);}
  unsigned long  operator++()                   {return increment(1);}
  unsigned long  count() const                  {return _count;}
  unsigned long  expCount() const               {return _expCount;}

private:
  void displayTic();

  std::ostream &    _os;  // may not be present in all imps
  const std::string _s1;  // string is more general, safer than 
  const std::string _s2;  // const char *, and efficiency or size are
  const std::string _s3;  // not issues

  unsigned long _count, _expCount, _nextTicCount;
  unsigned int  _tic;
};

#endif
