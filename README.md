imolecule
=========

An interactive 3D molecule viewer for the IPython notebook and
[browsers](http://www.patrick-fuller.com/imolecule/).

![](http://www.patrick-fuller.com/wp-content/uploads/2013/03/nu_100_1920.png)

Usage
=====

###Dependency

This requires the [Open Babel](http://openbabel.org/wiki/Main_Page) library
and Python bindings for chemical file format parsing. which can be installed
via apt-get/macports and pip.

```bash
apt-get/port install openbabel
pip install openbabel-python
```

In the case of OSX + homebrew, the default installer doesn't include python
bindings. Instead, use

```bash
brew install https://raw.github.com/rwest/homebrew/open-babel-new/Library/Formula/eigen2.rb
brew install https://raw.github.com/rwest/homebrew/open-babel-new/Library/Formula/open-babel.rb
```

If you have problems, refer to the Open Babel website.

###IPython notebook

Open a new notebook with `ipython notebook` and make sure that the `imolecule`
directory is either in the directory you started the notebook or your
PYTHONPATH. You can test the setup by typing:

```python
import imolecule
imolecule.draw("c1ccccc1")
```

into a notebook cell. You should get the output:

![](http://www.patrick-fuller.com/wp-content/uploads/2013/03/benzene_300.png)

The drawer can handle any format specified [here](http://openbabel.org/docs/2.3.1/FileFormats/Overview.html),
and can be set up to better handle different use cases. For example, to view a
crystallographic information file with pre-generated coordinates, you can
disable location generation and increase the size of the visualizer.

![](http://www.patrick-fuller.com/wp-content/uploads/2013/03/nu125_600.png)

###Full Browser

A version of the browser can be found at http://www.patrick-fuller.com/imolecule/.
To start your own local version, cd to the `imolecule` directory and start a
server with:

```bash
python -m SimpleHTTPServer
```

Navigate a browser to http://localhost:8000/index.html, and you're done.

The site allows for loading molecules via a simple file drag-and-drop interface.
The input takes json molecule files, which can be generated by using the
`molecule_to_json.py` script as a command-line tool:

```bash
python molecule_to_json.py smi c1ccccc1 > benzene.json
python molecule_to_json.py cif ../NU-125.cif --no-optimize > nu_125.json
```

Note that this requires the Open Babel library, as stated above.

