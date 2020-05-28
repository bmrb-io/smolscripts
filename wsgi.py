#!/usr/bin/env python
#
# need env in conda setup to get conda's python
# this is written for python 3.6.5
#

import os
import sys
import datetime
import pprint
import io

import werkzeug.routing
import werkzeug.wrappers
import werkzeug.exceptions
import werkzeug.utils

import bmrbmb.chemcomp

#
#
class App( object ) :
    #
    #
    def __init__( self ) :
        self._map = werkzeug.routing.Map ( [
            werkzeug.routing.Rule( "/help", endpoint = "help" ),
            werkzeug.routing.Rule( "/", endpoint = "new" ),
            werkzeug.routing.Rule( "/crunch", endpoint = "crunch" )
        ] )

        lastmod = datetime.datetime.fromtimestamp( os.path.getmtime( __file__ ) )
        lastmod = lastmod.replace( microsecond = 0 )
        self._lastmod = lastmod.isoformat( " " )

        self._mol = None

    ###################################################
    # werkzeug wsgi scaffolding
    #
    def __call__( self, environ, start_response ) :
        return self.wsgi_app( environ, start_response )

    #
    def wsgi_app( self, environ, start_response ) :
        request = werkzeug.wrappers.Request( environ )
        response = self.dispatch_request( request )
        return response( environ, start_response )

    #
    def dispatch_request( self, request ) :
        adapter = self._map.bind_to_environ( request.environ )
        try :
            endpoint, values = adapter.match()
            return getattr( self, "on_" + endpoint )( request, **values )
        except werkzeug.exceptions.HTTPException as e :
            return e

    #################################################
    # URL handlers
    #
    #
    def on_help( self, request ) :
        s = [ "<!doctype html>", '<html lang="en_us">', "<head>", '<meta charset="UTF-8">',
            "<title>molecule pix</title>", "</head>",
            "<body>\n",
            '<div style="padding: 2em; border: solid 1px black; border-radius: 3px; display: inline-block; font-size: 1.5em;">\n',
            "<ol>", "<li>Upload a MOL/SDF file,</li>", "<li>press a button,</li>", "<li>Profit.</li>", "</ol>\n",
            "</body></html>\n"]

        return werkzeug.wrappers.Response( s, status = 200, content_type = "text/html" )

    #
    #
    def on_new( self, request ) :

        s = [ "<!doctype html>", '<html lang="en_us">', "<head>", '<meta charset="UTF-8">',
            "<title>molecule pix</title>", "</head>",
            "<body>\n",
            '<div style="padding: 2em; border: solid 1px black; border-radius: 3px; display: inline-block; font-size: 1.25em;">\n',
            '<p style="text-align: center;">Crunch an SDF</p><hr size="1">\n' ]
        if not request.environ["SCRIPT_NAME"] == "" :
            if request.environ["SCRIPT_NAME"].endswith( "/" ) :
                s.append( '<form action="' + request.environ["SCRIPT_NAME"] + 'crunch" method="post" enctype="multipart/form-data">\n' )
            else :
                s.append( '<form action="' + request.environ["SCRIPT_NAME"] + '/crunch" method="post" enctype="multipart/form-data">\n' )
        else :
            s.append( '<form action="crunch" method="post" enctype="multipart/form-data">\n' )

        s.extend( [ '<p>Upload a MOL/SDF file: <input type="file" name="file" size="25"></p>\n',

# buttons are "pix" (2) and  "list"
#
            '<p style="text-align: center">\n<input type="submit" name="pix" value="create SVG (OpenBabel)">\n',
            '<input type="submit" name="rdpix" value="create SVG (RDKit)">\n',
            '<input type="submit" name="list" value="create atom list"></p>\n',
            'Minimize: <input type="radio" name="minimize" id="minimize" value="UFF">UFF',
            '<input type="radio" name="minimize" id="minimize" value="MMFF94">MMFF94',
            '<input type="radio" name="minimize" id="minimize" value="MMFF94s">MMFF94s\n',
            '</form>',
            '<p>If all else fails, try it in <a href="https://marvinjs-demo.chemaxon.com/latest/demo.html">Marvin</a></p>\n',
            "</div></body></html>\n"] )
        return werkzeug.wrappers.Response( s, status = 200, content_type = "text/html" )

    #
    #
    def on_crunch( self, request ) :

#        pprint.pprint( request.environ, stream = sys.stderr )
#        sys.stderr.write( "**************\n" )
#        pprint.pprint( request.files.values(), stream = sys.stderr )

        buf = io.StringIO()
        try :
            fstream = request.files.values().__next__()
        except AttributeError :
            fstream = request.files.values()[0]
        for line in fstream.stream : buf.write( bytes( line ).decode( "ascii" ) )
        sdf = buf.getvalue()
        buf.close()

        todo = "obpix"
        if "rdpix" in request.values.keys() : todo = "rdpix"
        elif "list" in request.values.keys() : todo = "list"
        request.close()

        if len( sdf ) < 1 :
            return werkzeug.wrappers.Response( ["Upload a file first"], status = 503, content_type = "text/html" )

        mol = bmrbmb.chemcomp.Molecule.from_string( string = str( sdf ) )

        if todo == "obpix" :
            return werkzeug.wrappers.Response( [mol._obmol._to_svg_full()], status = 200, content_type = "image/svg+xml" )

        if todo == "rdpix" :
            if "minimize" in request.values.keys() : mini = request.values["minimize"]
            else : mini = None
            return werkzeug.wrappers.Response( [mol._rdmol._to_svg_full( newstyle = True, minimize = mini )], status = 200, content_type = "image/svg+xml" )

        if todo == "list" :
            atoms = []
            for i in mol.iter_atoms() :
                atoms.append( str( i["id"] ) + "\n" )
            return werkzeug.wrappers.Response( atoms, status = 200, content_type = "text/plain" )

        return werkzeug.wrappers.Response( ["<!doctype html>", '<html lang="en_us">', "<head>", '<meta charset="UTF-8">',
            "<title>what now</title>", "</head>", "<body>\n",
            '<div style="padding: 2em; border: solid 1px black; border-radius: 3px; display: inline-block; font-size: 1.25em;">\n',
            "What now?"], status = 503, content_type = "text/html" )

####################################################################################################

#
# wsgi entry point
#
application = App()

#
#
#
if __name__ == '__main__':

    from werkzeug.serving import run_simple
    from werkzeug.debug import DebuggedApplication
    app = App()
    run_simple( '127.0.0.1', 5000, app, use_debugger = True, use_reloader = True )
