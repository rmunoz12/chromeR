# chromeR
A simple R package that processes identity by descent (IBD) segment data, such as that output by [ersa](https://github.com/rmunoz12/ersa), into a data set containing the lengths (in basepairs) of shared and not shared regions on each human autosome. Two functions are provided:

- `chrome_map`, which processes the IBD segment data, and
- `chrome_plot`, which shows an example plot of the output from `chrome_map`.

More information is provided in the help pages for each function.

## Installation
```
devtools::install_github("rmunoz12/chromeR")
```
