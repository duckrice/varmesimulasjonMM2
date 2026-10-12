@echo off
rem Dobbeltklikk for å åpne varmesimuleringen
cd /d "%~dp0"
start "" ".venv\Scripts\pythonw.exe" simulering.py
