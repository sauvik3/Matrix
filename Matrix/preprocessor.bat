@echo on
@for %%X in (cl.exe) do @ (set CLFOUND=%%~$PATH:X)
@if not defined CLFOUND (
    call "C:\Program Files\Microsoft Visual Studio 14.0\VC\vcvarsall.bat"
)

@cl.exe /P /D "MATRIX_USE_DLL" /D "DLL_EXPORTS" /I ".\include" /c .\src\matrix.cpp