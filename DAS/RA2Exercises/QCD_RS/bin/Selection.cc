#include "Selection.h"
#include <iostream>


void Cut(std::vector<Event*>& evts, float(*var)(const Event*), const char op, const float val)
{
  unsigned n=evts.size();
  if (op!='<'&&op!='>') return;
  for (std::vector<Event*>::iterator it=evts.begin();it<evts.end();){
    if (op=='<') {
      if ( var(*it)<val)
        it = evts.erase( it );
      else
        ++it;	
    } else if (op=='>'){
      if ( var(*it)>val)
        it = evts.erase( it );
      else
        ++it;	
    }   
  }
  std::cout<<"Removed "<<n-evts.size()<<" events, using func_ptr=&"<<var<<" "<<op<<" "<<val<<std::endl;
}
