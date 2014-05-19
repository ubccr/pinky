"""
FROWNS LICENSE

Copyright (c) 2001-2003, Brian Kelley
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met: 

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer. 
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following disclaimer
      in the documentation and/or other materials  provided with the
      distribution. 
    * Neither the name of Brian Kelley nor the names of frowns
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
import UserList, sys

class TokenHandler:
    def begin(self):
        """Called before any other callbacks"""
        pass
    def add_token(self, name, pos, text):
        """Called for each token; state name, character position, token text"""
        pass
    def error(self, errmsg, pos, text):
        """Called when an parse error occurs; error position, unparsed text"""
        raise "%s at position %d: %s" % (errmsg, pos, repr(text))
    def end(self):
        pass

class WriteHandler(TokenHandler):
    def __init__(self, outfile = sys.stdout):
        self.outfile = outfile
    def add_token(self, name, pos, text):
        self.outfile.write(" %d --> %s %s\n" % (pos, name, text))

class SaveTokens(TokenHandler, UserList.UserList):
    def begin(self):
        self.data[:] = []
    def add_token(self, name, pos, text):
        self.data.append( (name, pos, text) )

# A mixin class
class SilentErrors:
    def error(self, errmsg, pos, text):
        pass
