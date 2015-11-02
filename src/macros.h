#ifdef DEBUG
#  define ENTERSEXITS !Just for now (will be -D'ed)
#endif

#ifdef ENTERSEXITS
#  define ENTERS(routinename,err,error,linenum) \
          CALL Enters(routinename,err,error,linenum)
#  define EXITS(routinename) \
	  CALL Exits(routinename)
#  define ERRORS(routinename,err,error) \
	  CALL Errors(routinename,err,error)
#  define ERRORSEXITS(routinename,err,error) \
	  CALL Errors(routinename,err,error); CALL Exits(routinename)
#else
#  define ENTERS(routinename,err,error,linenum) !Do nothing
#  define EXITS(routinename) !Do nothing
#  define ERRORS(routinename,err,error) \
	  CALL Errors(routinename,err,error)
#  define ERRORSEXITS(routinename,err,error) \
	  CALL Errors(routinename,err,error)
#endif
