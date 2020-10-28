# LoopDetect - comprehensive detection of feedback loops in ODE models

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
LoopDetect is tested for Python 3, especially with Python version 3.8, but could also run with older Python versions.

In addition, LoopDetect can be found on [GitLab](https://gitlab.com/kabaum/loopdetect).



## Workflow and documentation

Function documentation is available on the LoopDetect pages, https://kabaum.gitlab.io/loopdetect/. There, you can also find a detailed [workflow description](https://kabaum.gitlab.io/loopdetect/workflow.html).


## Licensing
All code is licensed under the 3-clause BSD license, LoopDetect, Copyright (C) 2020  Katharina Baum.
