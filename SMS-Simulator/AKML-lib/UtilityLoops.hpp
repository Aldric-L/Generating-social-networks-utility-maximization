//
//  UtilityLoops.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 14/10/2023.
//

#ifndef UtilityLoops_h
#define UtilityLoops_h
namespace akml {
    template <int from, int to>
    struct for_ {
        template<template<int, int> class Fn, typename typearg>
        static void run(typearg arg) {
            Fn<from, from+1>::run(arg);
            for_<from + 1, to>::template run<Fn, typearg>(arg);
        }
        
        template<template<int> class Fn, typename typearg>
        static void run(typearg arg) {
            Fn<from>::run(arg);
            for_<from + 1, to>::template run<Fn, typearg>(arg);
        }
        
        template<template<int, int> class Fn, typename... typeargs>
        static void run(typeargs... args) {
            Fn<from, from+1>::run(args...);
            for_<from + 1, to>::template run<Fn, typeargs...>(args...);
        }
        
        template<template<int> class Fn, typename... typeargs>
        static void run(typeargs... args) {
            Fn<from>::run(args...);
            for_<from + 1, to>::template run<Fn, typeargs...>(args...);
        }
    };

    template <int to>
    struct for_<to, to> {
        template<template<int, int> class Fn, typename... typeargs>
        static void run(typeargs... args) {}
        template<template<int> class Fn, typename... typeargs>
        static void run(typeargs... args) {}
        
        template<template<int, int> class Fn, typename typearg>
        static void run(typearg arg) {}
        template<template<int> class Fn, typename typearg>
        static void run(typearg arg) {}
    };

}


#endif /* UtilityLoops_h */
