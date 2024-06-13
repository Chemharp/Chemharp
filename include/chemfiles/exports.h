
#ifndef CHFL_EXPORT_H
#define CHFL_EXPORT_H

#ifdef CHFL_STATIC_DEFINE
#  define CHFL_EXPORT
#  define CHFL_NO_EXPORT
#else
#  ifndef CHFL_EXPORT
#    ifdef chemfiles_EXPORTS
        /* We are building this library */
#      define CHFL_EXPORT 
#    else
        /* We are using this library */
#      define CHFL_EXPORT 
#    endif
#  endif

#  ifndef CHFL_NO_EXPORT
#    define CHFL_NO_EXPORT 
#  endif
#endif

#ifndef CHFL_DEPRECATED
#  define CHFL_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef CHFL_DEPRECATED_EXPORT
#  define CHFL_DEPRECATED_EXPORT CHFL_EXPORT CHFL_DEPRECATED
#endif

#ifndef CHFL_DEPRECATED_NO_EXPORT
#  define CHFL_DEPRECATED_NO_EXPORT CHFL_NO_EXPORT CHFL_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef CHFL_NO_DEPRECATED
#    define CHFL_NO_DEPRECATED
#  endif
#endif

#endif /* CHFL_EXPORT_H */
