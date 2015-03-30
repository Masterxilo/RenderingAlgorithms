::@echo off
setlocal enabledelayedexpansion

:: Compile
del /F rt\Main.class
del /F rt\AllTests.class
del /F rt\UnitTests.class

echo ====================================================
javac -cp *;. rt/Main.java rt/AllTests.java rt/UnitTests.java
echo ====================================================
echo ERRORLEVEL %ERRORLEVEL%

if %ERRORLEVEL% NEQ 0 (
goto end
)

If Not Exist "rt\Main.class" (
echo rt\Main.class not created
goto end
)
If Not Exist "rt\AllTests.class" (
echo rt\AllTests.class not created
goto end
)
If Not Exist "rt\UnitTests.class" (
echo rt\UnitTests.class not created
goto end
)

::if ERRORLEVEL 0 (
::if "%ERRORLEVEL%"=="0" (
mkdir output

:: Main, create output

echo ====================================================
java -ea -cp *;. rt.Main
echo ====================================================
echo ERRORLEVEL %ERRORLEVEL%

:: Unit tests, verify output
echo ====================================================
java -ea -cp *;. org.junit.runner.JUnitCore rt.AllTests
echo ====================================================
echo ERRORLEVEL %ERRORLEVEL%
::)

:end
