@echo off
rem Dobbeltklikk for å åpne varmesimuleringen.
rem Første gang lages et eget Python-miljø (.venv) med pakkene simuleringen trenger.
chcp 65001 >nul
cd /d "%~dp0"

if exist ".venv\installert" goto start

echo Første oppstart: installerer pakkene simuleringen trenger. Det tar et par minutter.
echo.
if not exist ".venv\Scripts\python.exe" (
    py -m venv .venv 2>/dev/null || python -m venv .venv
)
if not exist ".venv\Scripts\python.exe" (
    echo.
    echo Fant ikke Python. Installer Python 3.10 eller nyere fra https://www.python.org og prøv igjen.
    pause
    exit /b 1
)
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 (
    echo.
    echo Installasjonen feilet. Se meldingene over.
    pause
    exit /b 1
)
echo ok> ".venv\installert"

:start
start "" ".venv\Scripts\pythonw.exe" simulering.py
