# Documentation framework

The current documentation framework is built with `JupyterBook`. See [JupyterBook Documentation](https://jupyterbook.org/en/stable/intro.html) for detailed reference.

## JupyterBook setup
{bdg-secondary}`stable | Jupyter Book 0.15.1`

Install jupyter-book via python:
```bash
sudo python3 -m pip install jupyter-book
#
jupyter-book --version
```

## Building the document locally

Clean and build html pages locally with:
```bash
cd <REPOSITORY_PATH>
jupyter-book clean doc; jupyter-book build doc
```
The pages are built in `<REPOSITORY_PATH>/doc/_build/html/main/intro.html`.

## Doxygen reports

With doxygen, one can create code-level documentation. Follow the commands below:
```bash
# installation
sudo apt-get install doxygen doxygen-doc doxygen-gui graphviz 
# run doxygen
cd <REPOSITORY_PATH>/doc/doxygen
doxygen Doxyfile
cd doxyout/latex
make
```

## Add a new document to JupyterBook?

Find a suitable location in the folder structure. Add your new document with name `example.md` for simple markdown or `example.ipynb` for interactive notebook. Connect your document to the book in `_toc.yml` according to the structure. **Los gehts!**