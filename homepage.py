import os
from Flask

@app.route("/")
def homepage():
  contacts_dict={'Github':, 'https://github.com/Montmorency',
                 'email':'henry.lambert[cat.remove(c)]kcl.ac.uk',
                 '@condensedhank':'https://twitter.com/condensedhank'}
  cv_file = os.path.join(app.config['DOWNLOAD_FOLDER'], 'henry_lambert_cv.pdf')
  return render_template("homepage.html", contacts=contacts, cv_file=cv_file)

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
                  'Plasmonic Polarons':'http://link.aps.org/doi/10.1103/PhysRevLett.114.146404',
                  'Electron-phonon interaction in CaC6':'http://www.nature.com/articles/srep21414',
                  'Combined GW and Cumulant Expansion':'http://dx.doi.org/10.1103/PhysRevB.94.035103'
                 }
  return render_template("publications.html", publications=publications)

@app.route("/notes/")
def notes():
  """
  View of any notes, ideally in pdf format, that you may wish to share with the people.
  """
  notes = os.listdir(app.config['DOWNLOAD_FOLDER'])
  return render_template("notes.html", notes=notes)

@app.route('/interesting_stuff/')
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
                          "Total Perspective Clarity":"Total Perspective Clarity",
                          'Money Is The Way':'http://moneyistheway.blogspot.co.uk'}

  return render_template("interesting_stuff.html")

if __name__ == "__main__":
  SECRET_KEY = 'development key'
  DEBUG      = True
  DOWNLOAD_FOLDER = './download'
  ALLOWED_EXTENSIONS = set(['pdf'])
  app.config.from_object(__name__)    
  app.run()
