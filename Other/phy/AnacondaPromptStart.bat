@echo off
Title Anaconda-Prompt bat
Color 0B
Call :Browse4Folder "Choose your data:" "c:\scripts"
echo You have chosen this location "%Location%"
cd /d %Location%
REM echo %cd%
REM if %Location% == "Dialog Cancelled" (goto:eof)
call C:\Users\%USERNAME%\AppData\Local\Continuum\Anaconda3\Scripts\activate.bat 
REM call C:\Users\%USERNAME%\AppData\Local\Continuum\Anaconda3
REM call activate phy2
REM cmd/k "phy template-gui params.py"
cmd/k
exit
::***************************************************************************
:Browse4Folder
set Location=
set vbs="%temp%\_.vbs"
set cmd="%temp%\_.cmd"
for %%f in (%vbs% %cmd%) do if exist %%f del %%f
for %%g in ("vbs cmd") do if defined %%g set %%g=
(
    echo set shell=WScript.CreateObject("Shell.Application"^) 
    echo set f=shell.BrowseForFolder(0,"%~1",0,"%~2"^) 
    echo if typename(f^)="Nothing" Then  
    echo wscript.echo "set Location=Dialog Cancelled" 
    echo WScript.Quit(1^)
    echo end if 
    echo set fs=f.Items(^):set fi=fs.Item(^) 
    echo p=fi.Path:wscript.echo "set Location=" ^& p
)>%vbs%
cscript //nologo %vbs% > %cmd%
for /f "delims=" %%a in (%cmd%) do %%a
for %%f in (%vbs% %cmd%) do if exist %%f del /f /q %%f
for %%g in ("vbs cmd") do if defined %%g set %%g=
goto :eof
::***************************************************************************