import random, pdb, time
#red True, black False
#red black tree
class RBTree():
    def __init__(self, firstValue):
        self.top=RBTreeNode(firstValue, False)
    
    def insert(self, value):
        if value[0]>self.top.getValue()[0]:
            nextNode=self.top.getRight()
            if not nextNode:
                newNode=RBTreeNode(value, True, self.top, True)
                self.top.setRight(newNode)
                return
        else:
            nextNode=self.top.getLeft()
            if not nextNode:
                newNode=RBTreeNode(value, True, self.top, False)
                self.top.setLeft(newNode)
                return
        while True:
            if value[0]>nextNode.getValue()[0]:
                nextNode2=nextNode.getRight()
                if nextNode2:
                    nextNode=nextNode2
                else:
                    newNode=RBTreeNode(value, True, nextNode, True)
                    nextNode.setRight(newNode)
                    break
            else:
                nextNode2=nextNode.getLeft()
                if nextNode2:
                    nextNode=nextNode2
                else:
                    newNode=RBTreeNode(value, True, nextNode, False)
                    nextNode.setLeft(newNode)
                    break
        self.repairTree(newNode)
        
    def toList(self):
        l=[]
        self.dfs(self.top, l)
        return l
    
    def dfs(self, node, l):
        left_child=node.getLeft()
        right_child=node.getRight()
        if left_child:
            self.dfs(left_child, l)
        l.append(node.getValue())
        if right_child:
            self.dfs(right_child, l)
        
    def repairTree(self, current):
        #pdb.set_trace()
        parent=current.getParent()
        if parent:
            #red parent
            if parent.getColor():
                grandparent=parent.getParent()
                if parent.getSide():
                    sibling=grandparent.getLeft()
                else:
                    sibling=grandparent.getRight()
                #red uncle (case 3)
                if sibling and sibling.getColor():
                    parent.flipColor()
                    sibling.flipColor()
                    grandparent.setColor(True)
                    self.repairTree(grandparent)
                #black uncle (case 4)
                else:
                    #Downwards: Right, Left
                    if (not current.getSide()) and parent.getSide():
                        self.rotateRight(parent, current, current.getRight())
                        temp=parent
                        parent=current
                        current=temp
                    #Downwards: Left, Right
                    elif (not parent.getSide()) and current.getSide():
                        self.rotateLeft(parent, current, current.getLeft())
                        temp=parent
                        parent=current
                        current=temp
                    if current.getSide():
                        self.rotateLeft(grandparent, parent, parent.getLeft())
                    else:
                        self.rotateRight(grandparent, parent, parent.getRight())
                    parent.setColor(False)
                    grandparent.setColor(True)
            #else case is black parent (case2), no action required
        #is top of tree (case 1)
        else:
            current.setColor(False)
    
    def rotateLeft(self, grandparent, parent, left_child):
        ggp=grandparent.getParent()
        if ggp:
            if grandparent.getSide():
                ggp.setRight(parent)
            else:
                ggp.setLeft(parent)
        else:
            self.top=parent
        parent.setSide(grandparent.getSide())
        parent.setParent(ggp)
        parent.setLeft(grandparent)
        grandparent.setSide(False)
        grandparent.setParent(parent)
        grandparent.setRight(left_child)
        if left_child:
            left_child.setParent(grandparent)
            left_child.setSide(True)
        
    def rotateRight(self, grandparent, parent, right_child):
        ggp=grandparent.getParent()
        if ggp:
            if grandparent.getSide():
                ggp.setRight(parent)
            else:
                ggp.setLeft(parent)
        else:
            self.top=parent
        parent.setSide(grandparent.getSide())
        parent.setParent(ggp)
        parent.setRight(grandparent)
        grandparent.setSide(True)
        grandparent.setParent(parent)
        grandparent.setLeft(right_child)
        if right_child:
            right_child.setParent(grandparent)
            right_child.setSide(False)
        
        
class RBTreeNode():
    def __init__(self, value, color, parent=None, side=False):
        self.value=value
        self.left=None
        self.right=None
        self.color=color
        self.parent=parent
        #False for left, true for right
        self.side=side
    
    def setLeft(self, left=None):
        self.left=left
        
    def setRight(self, right=None):
        self.right=right
    
    def getLeft(self):
        return self.left
    
    def getRight(self):
        return self.right
    
    def getValue(self):
        return self.value
    
    def getColor(self):
        return self.color
    
    def flipColor(self):
        self.color=not self.color
    
    def setColor(self, color):
        self.color=color
        
    def getParent(self):
        return self.parent
    
    def setParent(self, parent=None):
        self.parent=parent
    
    def getSide(self):
        return self.side
    
    def setSide(self, side):
        self.side=side
        
if __name__=="__main__":
    if 1==1:
        size=10000000
        numbers=[]
        for i in range(size):
            numbers.append(random.random())
        first=numbers[-1]
        start=time.time()
        tree=RBTree([first,0])
        for i in range(size-1):
            tree.insert([numbers[i],0])
        a=tree.toList()
        end=time.time()
        print("Tree",end-start)
        start=time.time()
        l=[]
        for i in range(size):
            l.append([numbers[i],0])
        l.sort()
        end=time.time()
        print("List",end-start)
        del tree
        del l