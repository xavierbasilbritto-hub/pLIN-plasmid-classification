@echo off
REM ============================================================
REM  pLIN Launcher for Windows
REM  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
REM  Double-click this file to install dependencies and launch pLIN
REM ============================================================
title pLIN - Plasmid Life Identification Number System
echo.
echo  ============================================
echo   pLIN: Plasmid Life Identification Number
echo   Hierarchical Plasmid Classification System
echo  ============================================
echo.

REM Check for Python
where python >nul 2>nul
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Python is not installed or not in PATH.
    echo.
    echo Please install Python 3.10+ from https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during installation.
    echo.
    pause
    exit /b 1
)

REM Check Python version
for /f "tokens=2 delims= " %%v in ('python --version 2^>^&1') do set PYVER=%%v
echo [OK] Found Python %PYVER%

REM Set working directory to script location
cd /d "%~dp0"
echo [OK] Working directory: %CD%

REM Create virtual environment if it doesn't exist
if not exist ".venv\Scripts\activate.bat" (
    echo.
    echo [SETUP] Creating virtual environment...
    python -m venv .venv
    if %ERRORLEVEL% neq 0 (
        echo [ERROR] Failed to create virtual environment.
        pause
        exit /b 1
    )
    echo [OK] Virtual environment created.
)

REM Activate virtual environment
call .venv\Scripts\activate.bat
echo [OK] Virtual environment activated.

REM Install/upgrade dependencies
echo.
echo [SETUP] Installing dependencies (this may take a moment on first run)...
pip install --quiet --upgrade pip
pip install --quiet -r requirements.txt
if %ERRORLEVEL% neq 0 (
    echo [ERROR] Failed to install dependencies.
    echo Please check requirements.txt and try again.
    pause
    exit /b 1
)
echo [OK] All dependencies installed.

REM Check for AMRFinderPlus (optional)
where amrfinder >nul 2>nul
if %ERRORLEVEL% equ 0 (
    echo [OK] AMRFinderPlus detected — AMR analysis available.
) else (
    echo [INFO] AMRFinderPlus not found — pLIN will run without AMR analysis.
    echo        Install via: conda install -c bioconda -c conda-forge ncbi-amrfinderplus
)

REM Launch the app
echo.
echo ============================================
echo  Launching pLIN web application...
echo  The app will open in your default browser.
echo  Press Ctrl+C in this window to stop.
echo ============================================
echo.
streamlit run plin_app.py --server.headless false --browser.gatherUsageStats false

pause
