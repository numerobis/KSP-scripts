import itertools

###
### Adapted from http://stackoverflow.com/questions/280243/python-linked-list
###
### I removed a number of functions that I want to avoid ever using.

class LinkedList(object):
    """Immutable linked list class."""

    def __new__(cls, l=[]):
        if isinstance(l, LinkedList): return l # Immutable, so no copy needed.
        i = iter(l)
        try:
            head = i.next()
        except StopIteration:
            return cls.nil   # Return empty list singleton.

        tail = LinkedList(i)

        obj = super(LinkedList, cls).__new__(cls)
        obj._head = head
        obj._tail = tail
        return obj

    @classmethod
    def cons(cls, head, tail):
        ll =  cls([head])
        if not isinstance(tail, cls):
            tail = cls(tail)
        ll._tail = tail
        return ll

    # head and tail are not modifiable
    @property
    def head(self): return self._head

    @property
    def tail(self): return self._tail

    def __nonzero__(self): return True

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        x=self
        while x:
            yield x.head
            x=x.tail

    def __repr__(self):
        return "LinkedList([%s])" % ', '.join(map(repr,self))

class EmptyList(LinkedList):
    """A singleton representing an empty list."""
    def __new__(cls):
        return object.__new__(cls)

    def __iter__(self): return iter([])
    def __nonzero__(self): return False

    @property
    def head(self): raise IndexError("End of list")

    @property
    def tail(self): raise IndexError("End of list")

# Create EmptyList singleton
LinkedList.nil = EmptyList()
del EmptyList

nil = LinkedList.nil
cons = LinkedList.cons
