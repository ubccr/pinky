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


computes equivalence classes for atoms

equiv_class = compute_equiv_class(atom)
"""

def compute_equiv_class(atom):
    """(atom)->Computes a unique integer for an atom"""    
    try:
        equiv_class = atom.number + \
                      1000*(atom.charge+10) + \
                      100000*(atom.hcount) + \
                      1000000*(atom.weight)
    except TypeError:
        raise(ValueError, \
              "Can't compute number from atom.number %s atom.charge %s atom.hcount %s"\
              " atom.weight %s"%(atom.number, atom.charge, atom.hcount, atom.weight))
    return equiv_class
