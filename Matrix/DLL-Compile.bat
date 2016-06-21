@echo on
@for %%X in (cl.exe) do (set CLFOUND=%%~$PATH:X)
if not defined CLFOUND (
    call "C:\Program Files\Microsoft Visual Studio 14.0\VC\vcvarsall.bat"
)

@cl.exe /c /O2 /EHsc /MD /D "MATRIX_USE_DLL" /D "DLL_EXPORTS" /I ".\include" ".\src\matrix_exception.cpp" /Fo".\obj\matrix_exception.obj"
@cl.exe /c /O2 /EHsc /MD /D "MATRIX_USE_DLL" /D "DLL_EXPORTS" /I ".\include" ".\src\matrix.cpp" /Fo".\obj\matrix.obj"
@rc.exe -fo ".\rc\matrix.res" ".\rc\matrix.rc"

@link.exe /MACHINE:X86 /nodefaultlib:clib.lib /OPT:NOREF /DLL /subsystem:windows /out:".\bin\matrix.dll" ".\obj\matrix_exception.obj" ".\obj\matrix.obj" ".\rc\matrix.res"