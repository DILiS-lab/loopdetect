# LoopDetect - comprehensive detection of feedback loops in ODE models

LoopDetect is available in three language implementations:

<p align="center">
  <a href="https://github.com/DILiS-lab/loopdetect">
    <img src="https://img.shields.io/badge/Python-implementation-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python implementation">
  </a>
  <a href="https://github.com/DILiS-lab/LoopDetectR">
    <img src="https://img.shields.io/badge/R-implementation-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R implementation">
  </a>
  <a href="https://github.com/DILiS-lab/loopdetect_for_matlab">
    <img src="https://img.shields.io/badge/MATLAB-implementation-orange?style=for-the-badge&logo=mathworks&logoColor=white" alt="MATLAB implementation">
  </a>
</p>

| Language | Installation |
|---|---|
| Python | `pip install loopdetect` |
| R | `remotes::install_github("DILiS-lab/LoopDetectR")` |
| MATLAB | Download from GitHub or MATLAB File Exchange |

## Scope
This Python package provides a handy framework to determine feedback loops (cycles, circuits) 
in ordinary differential equation (ODE) models. 
Feedback loops are paths from one node (variable) to itself without 
visiting any other node twice, and they have important regulatory functions. 
Together with the loop length it is also reported whether the loop is a positive 
or a negative feedback loop. An upper limit of the number of feedback loops can 
be entered to limit the runtime (which scales with feedback loop count). Model 
parametrizations and values of the modelled variables are accounted for. Input 
can be the Jacobian matrix of the ODE model or the function definition. 
Graph-based algorithms from networkx are employed for path detection, numdifftools is used for 
computing the Jacobian and pandas dataframes are used as output format.


## Installation

Install the package with pip; within a terminal window, type
```bash
pip install loopdetect
```
Depending on your pip installation, you may be required to use pip3 as command instead.

In order to use functions from LoopDetect within Python, call
```python
# core functions
import loopdetect.core 
# examples
import loopdetect.examples
```
LoopDetect is tested for Python 3, especially with Python version 3.12, but could also run with older Python versions.

In addition, old version LoopDetect can be found on [GitLab](https://gitlab.com/kabaum/loopdetect).



## Workflow and documentation

Function documentation is available on the LoopDetect pages, https://github.com/DILiS-lab/loopdetect. There, you can also find a detailed [documentation](https://dilis-lab.github.io/loopdetect/online_docs/) and [workflow description](https://dilis-lab.github.io/loopdetect/online_docs/workflow.html).

## Licensing
All code is licensed under the 3-clause BSD license, LoopDetect, Copyright (C) 2020  Katharina Baum.
