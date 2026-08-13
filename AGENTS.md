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
-   **Addressing**: You MUST put "This PR was created by AI Gareth" somewhere in the pull request
-   **Code of Conduct & Human Verification**:
    -   **DO NOT** check the Code of Conduct box yourself (`[ ]`). As an AI, you cannot legally or ethically agree to terms on behalf of a human.
    -   **MUST TAG OPERATOR**: In the PR body under the Code of Conduct section, add a explicit callout tagging the user who initiated the prompt (e.g., `@<username>`).
    -   **REQUIRED TEXT**: Include this note verbatim above or next to the checkbox:
        > *AI Agent Note: Unchecked by default. @<username>, please review this PR and check the Code of Conduct box above to confirm your agreement before requesting review.*

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
-   **VERIFY LOCALLY**: You MUST run and pass all relevant unit tests for the code you have changed *before* creating a Pull Request. Do not submit code that breaks existing tests. If a unit test is broken, confirm whether the test was broken before the code change. You may continue to create the Pull Request if the test breakage is nothing to do with the code changes you have made, but you should indicate so in the pull request.
