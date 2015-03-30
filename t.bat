SETLOCAL ENABLEEXTENSIONS
set aaa=%~dpnx1
if "%aaa%"=="" (
	set aaa=all.java
)
javac %~dp0TangleHtml.java
javac %~dp0Tangle.java

if "%ERRORLEVEL%"=="0" (
java -ea -cp %~dp0 TangleHtml "%aaa%" %2 %3 %4 %5 %6 %7 %8 %9
java -ea -cp %~dp0 Tangle "%aaa%" %2 %3 %4 %5 %6 %7 %8 %9
goto end
) 
EXIT /B 1
:end
::debug
