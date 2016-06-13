def drepr(x, sort = True, indent = 0):
        """ Provides nice print format for a dictionnary """
        if isinstance(x, dict):
                r = '{\n'
                for (key, value) in (sorted(x.items()) if sort else x.iteritems()):
                        r += (' ' * (indent + 4)) + repr(key) + ': '
                        r += drepr(value, sort, indent + 4) + ',\n'
                r = r.rstrip(',\n') + '\n'
                r += (' ' * indent) + '}'
        # elif hasattr(x, '__iter__'):
        #       r = '[\n'
        #       for value in (sorted(x) if sort else x):
        #               r += (' ' * (indent + 4)) + drepr(value, sort, indent + 4) + ',\n'
        #       r = r.rstrip(',\n') + '\n'
        #       r += (' ' * indent) + ']'
        else:
                r = repr(x)
        return r
