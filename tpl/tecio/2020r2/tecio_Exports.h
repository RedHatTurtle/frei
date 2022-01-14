#ifndef TECPLOT_tecio_EXPORTS_H
#define TECPLOT_tecio_EXPORTS_H

/*
 * See http://gcc.gnu.org/wiki/Visibility for more information
 * about how this works.
 */

/* Generic helper definitions for shared library support */
#if defined _WIN32
    #define tecio_HELPER_DLL_IMPORT __declspec(dllimport)
    #define tecio_HELPER_DLL_EXPORT __declspec(dllexport)
    #define tecio_HELPER_DLL_LOCAL
#else
    #if __GNUC__ >= 4
        #define tecio_HELPER_DLL_IMPORT __attribute__ ((visibility("default")))
        #define tecio_HELPER_DLL_EXPORT __attribute__ ((visibility("default")))
        #define tecio_HELPER_DLL_LOCAL  __attribute__ ((visibility("hidden")))
    #else
        #define tecio_HELPER_DLL_IMPORT
        #define tecio_HELPER_DLL_EXPORT
        #define tecio_HELPER_DLL_LOCAL
    #endif
#endif

/*
 * The generic helper definitions above to define the following annotations:
 *   tecio_API:
 *       Used for functions and classes that should be part of the public API
 *       symbols. Methods and inner types of publicly declared API class
 *       symbols are also public API symbols.
 *   tecio_LOCAL:
 *       Used for methods and inner types of publicly declared API class
 *       symbols that should NOT be part of the public API symbols.
 * Any function or class symbol that is not annotated is considered private
 * unless it is a method or inner type of a publicly declared API class symbol.
 *
 * Note that static or template functions and classes should not be annotated.
 */
#ifdef tecio_EXPORTS // defined if we are building the tecio code (instead of using it)
  #define tecio_API tecio_HELPER_DLL_EXPORT
#else
  #define tecio_API tecio_HELPER_DLL_IMPORT
#endif // tecio_EXPORTS
#define tecio_LOCAL tecio_HELPER_DLL_LOCAL

/*
 * In VC++, if a class is flagged to be exported from a library, but is never used by that library,
 * it does not get exported. To force it to be exported, you must declare the class as exported
 * separately from the class declaration.
 *
 * >> ClassIWantExported.h
 *
 * class ClassIWantExported
 * {
 *   ...
 * };
 *
 * project_name_FORCE_EXPORT(ClassIWantExported);
 *
 */
#if defined _WIN32
    #define tecio_FORCE_EXPORT(ClassName) class tecio_API ClassName
#else
    #define tecio_FORCE_EXPORT(ClassName)
#endif

/*
 * In VC++, if two projects have classes that derive from a common templated base class, and
 * at least one of the projects exports its class publicly, it will cause multiple definition
 * link errors. The solution is to declare the base class as being exported as well.
 * http://msdn.microsoft.com/en-us/library/twa2aw10.aspx
 * LNK2005 and LNK1169 are commonly associated with this problem.
 *
 * >> Foo.h
 *
 * EXPORT_TEMPLATE_BASE(Bar<int>);
 *
 * class project_name_API Foo : public Bar<int>
 * {
 *   ...
 * };
 *
 */
#if defined _WIN32
    #define EXPORT_TEMPLATE_BASE(ClassName) template class __declspec(dllexport) ClassName;
#else
    #define EXPORT_TEMPLATE_BASE(ClassName)
#endif


#endif
