# AGENTS.md: System Context for AI Agents

This document provides essential context for AI agents assisting with development in the PyCBC repository. Adhere strictly to these instructions.

## 1. Top-Down Architecture Map

-   **Core Mission**: PyCBC is a software package for gravitational-wave astronomy, used to find and study binary merger signals in gravitational-wave data.
-   **Core Modules (`pycbc/`)**: The library is modular. Key components include:
    -   `pycbc.types`: Fundamental data structures (e.g., `TimeSeries`, `FrequencySeries`). This is the foundation for most data manipulation.
    -   `pycbc.psd`: Power Spectral Density (PSD) calculation, averaging, and manipulation.
    -   `pycbc.waveform`: Generation of gravitational waveform templates from physical models.
    -   `pycbc.filter`: Matched filtering, SNR calculation, and related signal processing algorithms.
    -   `pycbc.detector`: Detector geometry, antenna patterns, and time-of-flight calculations.
    -   `pycbc.inference`: A Bayesian inference framework for parameter estimation.
    -   `pycbc.workflow`: Pegasus workflow generation for orchestrating large-scale, distributed analyses on computing clusters.
-   **Entrypoints (`bin/`)**: Executable scripts that compose functionality from the `pycbc` library to perform analysis tasks (e.g., `pycbc_inspiral`, `pycbc_inference`). These are the primary user-facing tools.
-   **Core Primitives & Dependencies**:
    -   **LALSuite**: Heavy reliance on the LALSuite library for core gravitational-wave physics routines. **IMPORTANT**: There is a plan to remove this dependency. **LALSuite functions should be avoided as much as possible.** Note: There is a preference to use modern `igwn-ligolw` for I/O over older `lal` Python wrappers where possible.
    -   **GPU Acceleration**: Matched filtering can be accelerated using CUDA. See the `install` file for `PyCUDA` and related dependencies.
    -   **Data I/O**: Primary data format for banks and triggers is HDF5. Instrument data is in the GWF frame format. PSDs and configurations are often handled in XML or INI text files.

## 2. PR Creation & GitHub Workflow

-   **PR Labeling**: You MUST apply the `agent-assisted` label to any Pull Request you create or modify.
-   **PR Template**: You MUST locate and strictly adhere to the `.github/PULL_REQUEST_TEMPLATE.md` file. Populate all relevant sections of the template.
-   **Code of Conduct**: **CRITICAL**: Leave the Code of Conduct checkbox UNCHECKED: `[ ] The author of this pull request confirms they will adhere to the code of conduct`. As an AI, you cannot legally or ethically agree to this on behalf of a user.

## 3. Build, Test, and Quality Gates

-   **Environment Setup**:
    -   A LALSuite installation is a prerequisite. See the `install` file for details.
    -   Install PyCBC in editable mode from the root of the repository: `pip install -e .`
-   **Testing**:
    -   Run tests for the specific module you modified. For example, if you changed `pycbc/psd/read.py`, run: `pytest test/test_psd.py`.
    -   Before creating a PR, ensure all relevant tests pass locally.
-   **Code Quality**:
    -   This repository uses `flake8` for linting. The CI pipeline (`.github/workflows/check_code.yml`) specifically checks for unused imports (`F401`). Run `flake8 pycbc/ test/` to check your changes.
    -   The CI also runs `qlty` for static analysis. Ensure your code is clean to pass this gate.

## 4. Steering & Boundary Constraints (Hard Rules)

-   **DO NOT MODIFY**: Avoid modifying low-level C extensions, CUDA kernels, or files under `pycbc/lib` unless explicitly instructed by the user. These are performance-critical and require domain expertise.
-   **VERIFY LOCALLY**: You MUST run and pass all relevant unit tests for the code you have changed *before* creating a Pull Request. Do not submit code that breaks existing tests.
