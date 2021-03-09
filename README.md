# Web-server for Frag1D

## Description:
    This is the web-server implementation of Frag1D

    Frag1D is a method for predicting protein backbone diherial angles in 8 or
    3 states, given amino acid sequences of proteins.

    The web-server is developed with Django (>=2.0) and Python 3

    This software is open source and licensed under the MIT license (a copy
    of the license is included in the repository)


##Author
Nanjiang Shu

National Bioinformatics Infrastructure Sweden

Email: nanjiang.shu@scilifelab.se

## Reference

Tuping Zhou*, Nanjiang Shu* and Sven Hovmöller. A Novel Method for Accurate
One-dimensional Protein Structure Prediction based on Fragment Matching,
Bioinformatics, 2010;26(4):470-477. (*Co-first author)

## Installation

1. Install dependencies for the web server
    * Apache
    * mod\_wsgi

2. Install the virtual environments by 

    $ bash setup_virtualenv.sh

3. Create the django database db.sqlite3

4. Run 

    $ bash init.sh

    to initialize the working folder

5. In the folder `proj`, create a softlink of the setting script.

    For development version

        $ ln -s dev_settings.py settings.py

    For release version

        $ ln -s pro_settings.py settings.py

    Note: for the release version, you need to create a file with secret key
    and stored at `/etc/django_pro_secret_key.txt`

6.  On the computational node. run 

    $ virtualenv env --system-site-packages

    to make sure that python can use all other system-wide installed packages

