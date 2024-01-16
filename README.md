# KST toolbox

Toolbox authors:  Andrea Brancaccio, Luca Stefanutti &  Debora de Chiusole

Other toolbox references: Brancaccio A., de Chiusole D., Wickelmaier  F. (2024). Software. In Heller J. & Stefanutti L. *Knowledge structures Recent developments in Theory and Application.* Vol. 7. World Scientific. [DOI](https://doi.org/10.1142/13519)


The first MATLAB toolbox entirely devoted to knowledge structure theory.
Look at `GettingStarted.mlx` for more information on how to use the Toolbox.

### Data Structure

| Function  | Description                                                                             |
|------------------------------------|------------------------------------|
| `states`  | is a binary matrix $$k \times n$$ where $'k'$ is number of knowledge states and $'n'$ the number of problems in the domain $'Q'$. Each row $'0 < i \le k'$ is a vector representation of a knowledge state where a cell $'(i,j)= 1'$ if problem $'j'$ belongs to the knowledge states, 0 otherwise.
| `data` | is a binary matrix $'s \times n'$ where $'s'$ is number number of participant and $'n'$ the number of problems in the domain $'Q'$. Each row $'0 < i \le s'$ is a vector representation of participants $'i'$ response pattern. A cell $'(i,j)= 1'$ means that participant $'i'$ responded correctly at problem $'j'$, 0 otherwise. |
| `model` | is the main output from all functions that applied a version of the basic local Independence model. It is a structure with multiple fields corresponding to all the relevant features of the fitted models such as the parameters, absolute fit statistics, and model comparison statistics.|


### Building and plotting knowledge structures

| Function  | Description                                                                             |
|------------------------------------|------------------------------------|
| `skillmap`  |Delineate the knowledge structure from a skill map|
| `base2struct` | Build a knowledge space from its basis (algorithm described in Doignon & Falmagne, 1999, p. 32)|
| `plot_knowledge_structure` | Visualize the Hasse diagram of a knowledge structure |

### Fitting the Basic Local Indipendence Model

| Function  | Description                                                                             |
|------------------------------------|------------------------------------|
| `blim`  |Estimate the BLIM by maximum-likelihood via the EM algorithm|
| `sample` | Sample based on the BLIM’s parameters|
| `bootblim` | Compute the Chi-square p-value of the BLIM via parametric bootstrap |
| `blimfit` | Estimate the BLIM and compute its p-value via bootstrap|

### Testing the Local Identifiability

 Function  | Description             |
|------------------------------------|------------------------------------|
| `blimit`  |basic local independence model identification test see Stefanutti, L., Heller, J., Anselmi, P., and Robusto, E. (2012).  |

<p xmlns:cc="http://creativecommons.org/ns#" xmlns:dct="http://purl.org/dc/terms/"><span property="dct:title">KST-toolbox</span> is licensed under <a href="http://creativecommons.org/licenses/by-nc-sa/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY-NC-SA 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/nc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/sa.svg?ref=chooser-v1"></a></p>