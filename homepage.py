import os
import json
from flask import Flask, render_template, request, url_for, flash, g, session, send_file

app = Flask(__name__)

@app.route("/")
def homepage():
  contacts={'Github':'https://github.com/Montmorency',
            'email':'lamberh@tcd.ie',
            'twitter':'https://twitter.com/condensedhank'}
  return render_template("homepage.html", contacts=contacts)

@app.route("/_download/<path:filename>")
def download_file(filename):
  file_path=os.path.join(app.config['DOWNLOAD_FOLDER'], filename) 
  return send_file(file_path)

@app.route("/publications/")
def publications():
  """
  View with links to published papers.
  """
  publications = {'Ab initio Sternheimer-GW method': 'http://link.aps.org/doi/10.1103/PhysRevB.88.075117',
                  'Plasmonic polarons':'http://link.aps.org/doi/10.1103/PhysRevLett.114.146404',
                  'Electron-phonon interaction in CaC6':'http://www.nature.com/articles/srep21414',
                  'Combined GW and cumulant expansion':'http://dx.doi.org/10.1103/PhysRevB.94.035103'}
  hems_publications = {'Hydrogen embrittlement: future directions: discussion':'http://rsta.royalsocietypublishing.org/content/375/2098/20170029',
                       'Hydrogen transport and trapping: from quantum effects to alloy design: discussion':'http://rsta.royalsocietypublishing.org/content/375/2098/20170031'}
  return render_template("publications.html", publications=publications, hems_publications=hems_publications)

@app.route("/sorb/")
def sorb():
  from sorb.sorb import sorb_extracts
  print sorb_extracts
  return render_template("sorb.html", sorb_extracts=sorb_extracts)

@app.route("/environment/")
def environment():
  return render_template("environment.html")

@app.route("/notes/")
def notes():
  """
  View of any notes, ideally in pdf format, that you may wish to share with the people.
  """
  notes = os.listdir(app.config['DOWNLOAD_FOLDER'])
  return render_template("notes.html", notes=notes)

@app.route("/data_structures/")
def data_structures():
  """
  Visualizations of different sorting algorithms.
  """
  return render_template("data_structures.html")

@app.route("/quaternions")
def quaternions():
  """
  View for handling quaternion notes.
  """
  return render_template("quaternions.html")

@app.route('/tensorflow_notes/')
def tensorflow_notes():
  """
  WebNotes for using tensorflow.
  """
  return render_template('tensorflow_notes.html')

@app.route('/interesting_stuff/')
def interesting_stuff():
  """
  View of links to alternative pages.
  """
  interesting_websites = {'Tensorflow':'https://www.tensorflow.org',
                          'Quantum Espresso':'http://www.quantum-espresso.org',
                          'QUIP': "http://libatoms.github.io/QUIP/index.html",
                          "Mr. Bauld's English":'http://www.mrbauld.com',
                          "Flask Documentation":"http://flask.pocoo.org/docs/0.11/",
                          "Peewee ORM": "http://peewee.readthedocs.io/en/latest/",
                          'Learn You a Haskell For Great Good':'http://learnyouahaskell.com',
                          "Python":"https://www.python.org",
                          'Racing Post': 'http://www.racingpost.com',
                          "Sir David Mackay FRS":'http://www.inference.phy.cam.ac.uk/mackay/',
                          'Wolfson DVD Library':"https://www.wolfson.ox.ac.uk/dvd-library",
                          'Money Is The Way':'http://moneyistheway.blogspot.co.uk'}
  superconductivity = {"Philip Anderson: BCS Scientific Love of my Life.":"http://dx.doi.org/10.1142/S0217979210056426"}
  return render_template("interesting_stuff.html", interesting_websites=interesting_websites, superconductivity=superconductivity)

app.config['DEBUG']              = True
app.config['SECRET_KEY']         = '123412sdfalkjasflksqejvnoryyclzpiej'
app.config['DOWNLOAD_FOLDER']    = os.environ['DOWNLOAD_FOLDER']
app.config['ALLOWED_EXTENSIONS'] = set(['pdf'])

if __name__ == "__main__":
  app.run()
