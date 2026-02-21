@echo off
REM ═══════════════════════════════════════════════════════════════════════════════
REM  pLIN — All-in-One Setup & Launch for Windows
REM  Copyright (C) 2025 Basil Xavier Britto — GPL-3.0 + Citation clause
REM
REM  Single file: installs everything + launches the pLIN web application.
REM
REM  Usage:
REM    pLIN_Windows.bat                  Full install + launch (default)
REM    pLIN_Windows.bat --install        Install only (no launch)
REM    pLIN_Windows.bat --launch         Launch only (skip install)
REM    pLIN_Windows.bat --check          Check environment only
REM    pLIN_Windows.bat --uninstall      Remove conda env + venv
REM    pLIN_Windows.bat --docker         Build & run via Docker
REM
REM  Works on: Windows 10 (1903+), Windows 11
REM            x64 processors, WSL2 supported
REM ═══════════════════════════════════════════════════════════════════════════════

setlocal EnableDelayedExpansion
title pLIN v3.0.0 — Plasmid Lineage Identification Number System

REM ── Change to script directory ──────────────────────────────────────────────
cd /d "%~dp0"

REM ── Parse arguments ─────────────────────────────────────────────────────────
set "MODE=full"
if "%~1"=="--install"   set "MODE=install"
if "%~1"=="--launch"    set "MODE=launch"
if "%~1"=="--check"     set "MODE=check"
if "%~1"=="--uninstall" set "MODE=uninstall"
if "%~1"=="--docker"    set "MODE=docker"
if "%~1"=="--help"      set "MODE=help"
if "%~1"=="-h"          set "MODE=help"

REM ── Route to mode ───────────────────────────────────────────────────────────
if "%MODE%"=="help"      goto :show_help
if "%MODE%"=="docker"    goto :docker_mode
if "%MODE%"=="uninstall" goto :uninstall_mode
if "%MODE%"=="check"     goto :check_mode
if "%MODE%"=="launch"    goto :launch_only
if "%MODE%"=="install"   goto :install_only
if "%MODE%"=="full"      goto :full_mode

:show_help
call :banner
echo.
echo   Usage: pLIN_Windows.bat [OPTION]
echo.
echo   Options:
echo     (no option)    Full install + launch (default)
echo     --install      Install only (no launch)
echo     --launch       Launch only (skip install)
echo     --check        Check environment status
echo     --docker       Build and run via Docker
echo     --uninstall    Remove conda env + venv
echo     --help         Show this help message
echo.
echo   Requirements:
echo     Windows 10 (1903+) or Windows 11
echo     Python 3.10+ or Miniconda/Anaconda
echo.
echo   For full bioinformatics tool support, use WSL2:
echo     wsl --install
echo     Then run pLIN_Linux.sh inside WSL
echo.
goto :eof

REM ═══════════════════════════════════════════════════════════════════════════════
REM  FUNCTIONS
REM ═══════════════════════════════════════════════════════════════════════════════

:banner
echo.
echo  ================================================================
echo    pLIN: Plasmid Lineage Identification Number System
echo    Version 3.0.0 — All-in-One Installer for Windows
echo    Hierarchical Plasmid Classification + AMR Surveillance
echo  ================================================================
echo.
goto :eof

:ok
echo   [OK]   %~1
goto :eof

:warn
echo   [WARN] %~1
goto :eof

:fail
echo   [FAIL] %~1
goto :eof

:info
echo   [INFO] %~1
goto :eof

REM ── Detect platform ─────────────────────────────────────────────────────────
:detect_platform
echo   Platform:     Windows
for /f "tokens=2 delims==" %%a in ('wmic os get Caption /value 2^>nul ^| findstr Caption') do (
    echo   OS:           %%a
)
echo   Architecture: %PROCESSOR_ARCHITECTURE%
echo   Working dir:  %CD%
goto :eof

REM ── Find Python ─────────────────────────────────────────────────────────────
:find_python
set "PYTHON_BIN="
set "PYTHON_VER="

REM Check python in PATH
for %%p in (python python3) do (
    where %%p >nul 2>&1
    if !errorlevel! equ 0 (
        for /f "tokens=2 delims= " %%v in ('%%p --version 2^>^&1') do (
            set "PYTHON_VER=%%v"
        )
        for /f "tokens=1,2 delims=." %%a in ("!PYTHON_VER!") do (
            if %%a geq 3 if %%b geq 10 (
                set "PYTHON_BIN=%%p"
                goto :found_python
            )
        )
    )
)

REM Check conda environments
for %%d in (
    "%USERPROFILE%\miniconda3\python.exe"
    "%USERPROFILE%\Miniconda3\python.exe"
    "%USERPROFILE%\anaconda3\python.exe"
    "%USERPROFILE%\Anaconda3\python.exe"
    "%LOCALAPPDATA%\miniconda3\python.exe"
    "%PROGRAMDATA%\miniconda3\python.exe"
) do (
    if exist %%d (
        for /f "tokens=2 delims= " %%v in ('%%d --version 2^>^&1') do (
            set "PYTHON_VER=%%v"
        )
        set "PYTHON_BIN=%%~d"
        goto :found_python
    )
)

call :fail "Python 3.10+ not found"
call :info "Download from: https://www.python.org/downloads/"
call :info "IMPORTANT: Check 'Add Python to PATH' during installation"
set "PYTHON_BIN="
goto :eof

:found_python
call :ok "Python !PYTHON_VER! (!PYTHON_BIN!)"
goto :eof

REM ── Find conda ──────────────────────────────────────────────────────────────
:find_conda
set "CONDA_BIN="

where conda >nul 2>&1
if %errorlevel% equ 0 (
    for /f "delims=" %%p in ('where conda') do set "CONDA_BIN=%%p"
    goto :found_conda
)

REM Search common Windows locations
for %%d in (
    "%USERPROFILE%\miniconda3\Scripts\conda.exe"
    "%USERPROFILE%\Miniconda3\Scripts\conda.exe"
    "%USERPROFILE%\anaconda3\Scripts\conda.exe"
    "%USERPROFILE%\Anaconda3\Scripts\conda.exe"
    "%LOCALAPPDATA%\miniconda3\Scripts\conda.exe"
    "%PROGRAMDATA%\miniconda3\Scripts\conda.exe"
    "%USERPROFILE%\miniforge3\Scripts\conda.exe"
    "%USERPROFILE%\mambaforge\Scripts\conda.exe"
) do (
    if exist %%d (
        set "CONDA_BIN=%%~d"
        goto :found_conda
    )
)

set "CONDA_BIN="
goto :eof

:found_conda
call :ok "Conda: !CONDA_BIN!"
goto :eof

REM ── Check required files ────────────────────────────────────────────────────
:check_files
set "FILES_OK=1"
for %%f in (
    "plin_app.py"
    "requirements.txt"
    "data\inc_classifier.npz"
    "data\inc_centroids.npz"
    "output\pLIN_assignments.tsv"
) do (
    if exist "%%~f" (
        call :ok "Found %%~f"
    ) else (
        call :fail "Missing: %%~f"
        set "FILES_OK=0"
    )
)
goto :eof

REM ── Check Python packages ──────────────────────────────────────────────────
:check_packages
set "PKG_MISSING=0"
for %%p in (streamlit numpy pandas scipy Bio sklearn matplotlib seaborn plotly pptx requests joblib) do (
    !PYTHON_BIN! -c "import %%p" >nul 2>&1
    if !errorlevel! equ 0 (
        call :ok "Package: %%p"
    ) else (
        call :warn "Missing: %%p"
        set /a PKG_MISSING+=1
    )
)
goto :eof

REM ── Check bioinformatics tools ──────────────────────────────────────────────
:check_biotools
set "BIO_FOUND=0"
set "BIO_TOTAL=9"
for %%t in (amrfinder mash fastANI minimap2 minced blastn prodigal mob_typer mlst) do (
    where %%t >nul 2>&1
    if !errorlevel! equ 0 (
        call :ok "%%t"
        set /a BIO_FOUND+=1
    ) else (
        REM Check conda env
        set "TOOL_FOUND=0"
        for %%d in (
            "%USERPROFILE%\miniconda3\envs\pLIN_tools\Scripts\%%t.exe"
            "%USERPROFILE%\miniconda3\envs\pLIN_tools\bin\%%t"
            "%USERPROFILE%\Miniconda3\envs\pLIN_tools\Scripts\%%t.exe"
            "%USERPROFILE%\anaconda3\envs\pLIN_tools\Scripts\%%t.exe"
        ) do (
            if exist %%d (
                call :ok "%%t (in conda env)"
                set /a BIO_FOUND+=1
                set "TOOL_FOUND=1"
            )
        )
        if !TOOL_FOUND! equ 0 (
            call :warn "%%t not found (optional)"
        )
    )
)
goto :eof

REM ── Install Python packages ─────────────────────────────────────────────────
:install_python_packages
echo.
echo [%~1/%~2] Installing Python dependencies...

!PYTHON_BIN! -m pip install --upgrade pip --quiet >nul 2>&1
!PYTHON_BIN! -m pip install --quiet -r requirements.txt
if !errorlevel! equ 0 (
    call :ok "All Python packages installed"
) else (
    call :fail "Failed to install Python packages"
    call :info "Try manually: pip install -r requirements.txt"
)
goto :eof

REM ── Install biotools via conda ──────────────────────────────────────────────
:install_biotools
echo.
echo [%~1/%~2] Installing bioinformatics tools via conda...

call :find_conda
if "!CONDA_BIN!"=="" (
    call :warn "Conda not found — skipping bioinformatics tools"
    call :info "Core pLIN features work without these tools"
    call :info "Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html"
    call :info ""
    call :info "For full bioinformatics tool support on Windows, consider:"
    call :info "  1. Install WSL2:  wsl --install"
    call :info "  2. Run pLIN_Linux.sh inside WSL"
    goto :eof
)

REM Create or reuse conda env
"!CONDA_BIN!" env list 2>nul | findstr /C:"pLIN_tools" >nul 2>&1
if !errorlevel! neq 0 (
    call :info "Creating conda environment 'pLIN_tools' with Python 3.11..."
    "!CONDA_BIN!" create -n pLIN_tools python=3.11 -y --quiet >nul 2>&1
    call :ok "Conda environment created"
) else (
    call :info "Conda environment 'pLIN_tools' already exists"
)

call conda activate pLIN_tools

REM Configure channels
conda config --add channels defaults 2>nul
conda config --add channels bioconda 2>nul
conda config --add channels conda-forge 2>nul

REM Update PYTHON_BIN to use conda env python
set "PYTHON_BIN=python"

REM Install tools
set "INSTALLED=0"
for %%t in (ncbi-amrfinderplus mash fastani minimap2 minced blast prodigal mlst) do (
    call :info "Installing %%t..."
    conda install -y -c bioconda -c conda-forge %%t --quiet >nul 2>&1
    if !errorlevel! equ 0 (
        call :ok "Installed %%t"
        set /a INSTALLED+=1
    ) else (
        call :warn "Could not install %%t (optional)"
    )
)

REM MOBsuite via pip
call :info "Installing MOBsuite..."
pip install mob_suite --quiet >nul 2>&1
if !errorlevel! equ 0 (
    call :ok "Installed MOBsuite"
    set /a INSTALLED+=1
) else (
    call :warn "Could not install MOBsuite (optional)"
)

REM Update AMRFinderPlus database
where amrfinder >nul 2>&1
if !errorlevel! equ 0 (
    call :info "Updating AMRFinderPlus database..."
    amrfinder --update >nul 2>&1
)

call :ok "Installed !INSTALLED!/9 bioinformatics tools"
goto :eof

REM ── Setup directories ───────────────────────────────────────────────────────
:setup_directories
if not exist "output\amrfinder" mkdir "output\amrfinder"
if not exist "output\crispr_evaluation" mkdir "output\crispr_evaluation"
if not exist "output\figures" mkdir "output\figures"
if not exist "output\manuscripts\docx" mkdir "output\manuscripts\docx"
goto :eof

REM ── Launch app ──────────────────────────────────────────────────────────────
:launch_app
echo.
echo [%~1/%~2] Launching pLIN...

if not exist "plin_app.py" (
    call :fail "plin_app.py not found in current directory"
    pause
    goto :eof
)

REM Try to activate conda env
call :find_conda
if not "!CONDA_BIN!"=="" (
    "!CONDA_BIN!" env list 2>nul | findstr /C:"pLIN_tools" >nul 2>&1
    if !errorlevel! equ 0 (
        call conda activate pLIN_tools 2>nul
    )
)

REM Find streamlit
set "STREAMLIT="
where streamlit >nul 2>&1
if !errorlevel! equ 0 (
    set "STREAMLIT=streamlit"
) else (
    REM Check conda env
    for %%d in (
        "%USERPROFILE%\miniconda3\envs\pLIN_tools\Scripts\streamlit.exe"
        "%USERPROFILE%\Miniconda3\envs\pLIN_tools\Scripts\streamlit.exe"
        "%USERPROFILE%\anaconda3\envs\pLIN_tools\Scripts\streamlit.exe"
    ) do (
        if exist %%d (
            set "STREAMLIT=%%~d"
        )
    )
    REM Check venv
    if exist ".venv\Scripts\streamlit.exe" (
        set "STREAMLIT=.venv\Scripts\streamlit.exe"
    )
)

if "!STREAMLIT!"=="" (
    call :fail "Streamlit not found. Run: pLIN_Windows.bat --install"
    pause
    goto :eof
)

echo.
echo  ================================================================
echo    pLIN is starting at: http://localhost:8501
echo    Press Ctrl+C to stop the server.
echo  ================================================================
echo.

REM Open browser
start "" "http://localhost:8501" 2>nul

"!STREAMLIT!" run plin_app.py --server.headless true --server.port 8501 --server.address 0.0.0.0 --browser.gatherUsageStats false --server.maxUploadSize 200
goto :eof

REM ═══════════════════════════════════════════════════════════════════════════════
REM  MODES
REM ═══════════════════════════════════════════════════════════════════════════════

:check_mode
call :banner
call :detect_platform
echo.

echo.
echo [1/4] Checking Python...
call :find_python

echo.
echo [2/4] Checking required files...
call :check_files

echo.
echo [3/4] Checking Python packages...
if not "!PYTHON_BIN!"=="" (
    call :check_packages
) else (
    call :warn "Cannot check packages — Python not found"
)

echo.
echo [4/4] Checking bioinformatics tools...
call :check_biotools

echo.
echo  ================================================================
echo    Environment Summary
echo  ================================================================
echo    Bioinformatics tools: !BIO_FOUND!/!BIO_TOTAL! available
echo  ================================================================
echo.
pause
goto :eof

:docker_mode
call :banner

echo.
echo [1/3] Checking Docker...
where docker >nul 2>&1
if !errorlevel! neq 0 (
    call :fail "Docker not found"
    call :info "Install Docker Desktop: https://www.docker.com/products/docker-desktop/"
    pause
    goto :eof
)
call :ok "Docker found"

echo.
echo [2/3] Building pLIN Docker image...
docker build -t plin:latest .
if !errorlevel! neq 0 (
    call :fail "Docker build failed"
    pause
    goto :eof
)
call :ok "Docker image built: plin:latest"

echo.
echo [3/3] Starting pLIN container...
docker stop plin-app 2>nul
docker rm plin-app 2>nul
docker run -d --name plin-app -p 8501:8501 -v "%CD%\data:/app/data" -v "%CD%\output:/app/output" plin:latest
if !errorlevel! equ 0 (
    call :ok "pLIN running in Docker"
    echo.
    echo   Open in browser: http://localhost:8501
    echo.
    call :info "Stop with:  docker stop plin-app"
    call :info "Logs:       docker logs -f plin-app"
    start "" "http://localhost:8501"
) else (
    call :fail "Failed to start container"
)
echo.
pause
goto :eof

:uninstall_mode
call :banner
echo   This will remove the pLIN conda environment and virtual environment.
echo   Your data files and sequences will NOT be deleted.
echo.
set /p "CONFIRM=  Continue? [y/N]: "
if /i not "!CONFIRM!"=="y" (
    echo   Cancelled.
    goto :eof
)

call :find_conda
if not "!CONDA_BIN!"=="" (
    "!CONDA_BIN!" env list 2>nul | findstr /C:"pLIN_tools" >nul 2>&1
    if !errorlevel! equ 0 (
        call :info "Removing conda environment 'pLIN_tools'..."
        "!CONDA_BIN!" env remove -n pLIN_tools -y >nul 2>&1
        call :ok "Conda environment removed"
    )
)

if exist ".venv" (
    call :info "Removing virtual environment..."
    rmdir /s /q ".venv"
    call :ok "Virtual environment removed"
)

call :ok "Uninstall complete. Data files are preserved."
echo.
pause
goto :eof

:install_only
call :banner

echo.
echo [1/5] Checking environment...
call :detect_platform
echo.
call :find_python
if "!PYTHON_BIN!"=="" (
    call :warn "Python not found — will install via conda"
)
call :check_files
if "!FILES_OK!"=="0" (
    call :fail "Required data files are missing."
    call :fail "Ensure you have the complete pLIN distribution."
    pause
    goto :eof
)

call :install_biotools 2 5
call :install_python_packages 3 5

echo.
echo [4/5] Setting up directories...
call :setup_directories
call :ok "Output directories ready"

echo.
echo  ================================================================
echo    Installation Complete!
echo  ================================================================
echo.
echo   Launch pLIN with: pLIN_Windows.bat --launch
echo   Or directly:      conda activate pLIN_tools ^& streamlit run plin_app.py
echo.
echo   For full bioinformatics tool support, consider using WSL2:
echo     wsl --install
echo     Then run pLIN_Linux.sh inside WSL
echo.
pause
goto :eof

:launch_only
call :banner
call :find_python
call :launch_app 1 1
goto :eof

:full_mode
call :banner

echo.
echo [1/5] Checking environment...
call :detect_platform
echo.
call :find_python
if "!PYTHON_BIN!"=="" (
    call :warn "Python not found — will install via conda"
)
call :check_files
if "!FILES_OK!"=="0" (
    call :fail "Required data files are missing."
    call :fail "Ensure you have the complete pLIN distribution."
    pause
    goto :eof
)

call :install_biotools 2 5
call :install_python_packages 3 5

echo.
echo [4/5] Setting up directories...
call :setup_directories
call :ok "Output directories ready"

call :launch_app 5 5
goto :eof
