

def reroot(tree, lab):
    for node in tree:
        if node.label == lab:
            #print "Found" , lab
            tree = node
            tree.parent=None
            return tree
    raise lab + " Not found"
        
if __name__=="""__main__""":

    import timeit

    import phylotree
    import newick
    import branch_lengths



    tree = newick.read_trees(open("ecol-let-sim/globalphylo").read())[0]
    #tree = newick.read_trees(open("examples/ca-phylo.new").read())[0]
    ages = branch_lengths.get_age_dict("ecol-let-sim/ages.txt")


    branch_lengths.bl_one(tree)

    ntree = tree.copy()

    # print tree
    branch_lengths.bl_bladj(tree, ages)

    #print ntree.label
    #print ages[ntree.label]
    branch_lengths.bl_bladj(ntree, ages)


    #print len(tree.leaves())
    #tree = reroot(tree,"moraceae")
    #print tree.label, len(tree.leaves())
    #ntree = reroot(ntree,"moraceae")

    tree.normalize()
    ntree.normalize()

    #print ntree.write(True) + ";"
    
    old_nodeages =  tree.node_ages()
    new_node_ages=ntree.node_ages()

    assert len(old_nodeages) == len(new_node_ages)
            
    # for i in range(len(old_nodeages)):
    #     if abs(old_nodeages[i][1] - new_node_ages[i][1] ) > 0.001 :
    #         print old_nodeages[i][0],  old_nodeages[i][1] - new_node_ages[i][1]
        
    
    tbl1 = timeit.Timer("branch_lengths.bl_bladj(tree, ages)", "from __main__ import tree, ages; import branch_lengths")    
#   tbl2 = timeit.Timer("branch_lengths.bl_bladj2(tree, ages)", "from __main__ import tree, ages; import branch_lengths")

    print "bladj1", tbl1.timeit(5)
#    print "bladj2", tbl2.timeit(10)
    
    #print tree
    
    t= timeit.Timer("tree.write(True)","from __main__ import tree")
    t2= timeit.Timer("tree.write2(True)","from __main__ import tree")
    t3= timeit.Timer("tree.write3(True)","from __main__ import tree")
    t4= timeit.Timer("tree.write4(True)","from __main__ import tree")

    #print tree.write4(True)
    

#    print "write4",t4.timeit(1000)    
#    print "write", t.timeit(1000)
#    print "write3",t3.timeit(1000)
    
#    print "write2", min(t2.repeat(10,100))
#    print "write", min(t.repeat(10,100))    
#    print "write3", min(t3.repeat(10,100))
#    print "write4", min(t4.repeat(10,100))

    
    
