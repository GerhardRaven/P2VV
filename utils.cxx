#include "utils.h"
char *_Format_( const char *fmt, ... ) {
        va_list ap;
        int size = 128;
        char *p(0);
        do { 
            p = new char[size];
            if ( p == 0) return 0; // give up, no memory...
            va_start(ap,fmt);
            int s = vsnprintf(p,size,fmt,ap);
            va_end(ap);
            if (s<size) break;
            // OOPS -- we truncated... delete p, double size, try again...
            size *= 2;
            delete[] p;
            p = 0;
        } while ( size < 2048 ) ; // don't continue doubling forever
        return p;
}

