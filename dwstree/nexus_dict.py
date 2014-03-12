######################################################
# NexusDict Class
######################################################

class NexusDict(dict):
    """Dictionary, that has case-insensitive keys.
       Dictionary also maintains order of items.
    
    Keys are retained in their original form
    when queried with keys() or items().

    Implementation: Inherites from dict. All key lookups are done
    against the uppercase keys, but all methods that expose keys to
    the user retrieve the original keys."""
    
    def __init__(self,*args,**kwargs):
        """Inherit from dict, store values as (key,value) pairs. Store
        uppercase keys in _fields"""
        dict.__init__(self,*args,**kwargs)
        self._fields = []

    def __getitem__(self, key):
        """Retrieve the value associated with 'key' (in any case)."""
        return dict.__getitem__(self,key.upper())[1]

    def __setitem__(self, key, value):
        """Associate 'value' with 'key'. If 'key' already exists, but
        in different case, it will be replaced."""
        k = key.upper()
        dict.__setitem__(self, k, (key, value))
        if not self.has_key(key):
            self._fields.append(k)

    def __delitem__(self,key) :
        key = key.upper()
        try:
            dict.__delitem__(self, key)
        except KeyError:
            pass
        try:
            self._fields.remove(key)
        except ValueError:
            pass
    
    def has_key(self, key):
        """Case insensitive test if 'key' exists."""
        k = key.upper()
        return k in self._fields

    def keys(self):
        """List of keys in their original case."""
        return [dict.__getitem__(self,key.upper())[0] for key in self._fields]

    def values(self):
        """List of values."""
        return [dict.__getitem__(self,key.upper())[1] for key in self._fields]
        

    def items(self):
        """List of (key,value) pairs."""
        return [dict.__getitem__(self,key.upper()) for key in self._fields]

    def get(self, key, default):
        """Retrieve value associated with 'key' or return default value return
            default"""
        return dict.get(self, key.upper(), default)[1]
         

    def setdefault(self, key, default):
        """If 'key' doesn't exists, associate it with the 'default' value.
        Return value associated with 'key'."""
        if not self.has_key(key):
            self[key] = default
        return self[key]

    def __repr__(self):
        """String representation of the dictionary."""
        items = ", ".join([("%r: %r" % (k,v)) for k,v in self.items()])
        return "{%s}" % items

    def __str__(self):
        """String representation of the dictionary."""
        return repr(self)


