import Page
from bintrees import RBTree


class Element:
    def __init__(self, page: Page):
        self.coefMap = RBTree()
        self.page = page

    @property
    def bigrade(self):
        return
