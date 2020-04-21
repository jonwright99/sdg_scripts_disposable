def check_neighbourhood(node,next_node,distance,end_size):
    l_min=(distance+end_size-1000) * .9
    l_max=(distance+end_size+1000) * 1.1


    node_found=0
    long_enough=0
    next_found=0
    for rid in lrr.node_reads[abs(node)]:
        t=lrr.read_threads[rid]

        node_fw=[x for x in t if x.node==node]
        node_bw=[x for x in t if x.node==-node]
        if len(node_fw)==1 and len(node_bw)==0:
            node_found+=1
            nf=node_fw[0]
            if lords.get_read_size(rid)>nf.end+l_max:
                long_enough+=1
                if [x for x in t if x.node==next_node and l_min<x.start-nf.end<l_max ]:
                    next_found+=1
        if len(node_bw)==1 and len(node_fw)==0:
            node_found+=1
            nb=node_bw[0]
            if nb.start<l_max:
                long_enough+=1
                if [x for x in t if x.node==-next_node and l_min<nb.start-x.end<l_max ]:
                    next_found+=1
    return (len(lrr.node_reads[abs(node)]),node_found,long_enough,next_found)

def review_thread(t):
    for i in range(len(t)):
        for j in range(len(t)):
            if i==j: print("  --  \t",end="")
            else:
                if i<j:
                    ns=check_neighbourhood(t[i].node,t[j].node,t[j].start-t[i].end,500)
                else:
                    ns=check_neighbourhood(-t[i].node,-t[j].node,t[i].start-t[j].end,500)
                s=100.0*ns[3]/ns[2]
                print("%3.2f\t"%(s),end="")
        print()

#find "truly connected" node pairs (cleanup read threads to only include large unique nodes, this could use some tapping)
def get_1to1_connections(NODE_SIZE=500,MAX_KCI=1.5,WIN_PERC=.75,WIN_MIN=10):


    hs_nodes=set([nv.node_id() for nv in ws.sdg.get_all_nodeviews() if nv.size()>=NODE_SIZE and nv.kci()<=MAX_KCI])
    for rid in range(len(lrr.read_threads)):
        lrr.read_threads[rid]=[x for x in lrr.read_threads[rid] if abs(x.node) in hs_nodes]
    mldg_specific=lrr.dg_from_threads(False)

    conn_ends=set()

    for nvf in mldg_specific.get_all_nodeviews():
        #print(nvf)
        for nv in [nvf, nvf.rc()]:
            #print(nv)
            try:
                mc=Counter([x.node().node_id() for x in nv.next()]).most_common(1)[0]
            except:
                continue
            #print ("\nnode %s\nmc=%s"%(str(nv),mc))
            #print(nv.next())
            if mc[1]>=WIN_MIN and mc[1]/len(nv.next())>=WIN_PERC:
                n=mc[0]
                mc2=Counter([x.node().node_id() for x in mldg_specific.get_nodeview(n).prev()]).most_common(1)[0]
                #print("mc2=%d"%mc2)
                if mc2[0]==nv.node_id() and mc2[1]>=WIN_MIN and mc2[1]/len(mldg_specific.get_nodeview(n).prev())>=WIN_PERC:
                    conn_ends.add((min(-nv.node_id(),n),max(-nv.node_id(),n)))

    print(len(conn_ends),"connections detected between",len(hs_nodes),"long unique nodes")
    from statistics import median
    conns=[]
    for c in conn_ends:
        n1=-c[0]
        n2=c[1]
        #rids=[x.support().id for x in mldg_specific.get_nodeview(n1).next() if x.node().node_id()==n2]
        dists=[x.distance() for x in mldg_specific.get_nodeview(n1).next() if x.node().node_id()==n2]
        conns.append([n1,n2,median(dists),len(dists)])

    return conns,mldg_specific

def solve_with_pf():
    ge=SDG.GraphEditor(ws)
    conn=get_1to1_connections()
    #for c in conn[0][2:3]:
    for c in conn[0]:
        print("\n\n",c)
        #for x in conn[1].get_nodeview(c[0]).next():
        #    print(x,x.support().id,len(lrr.read_threads[x.support().id]))
        if int(c[2])<-50 and [x for x in ws.sdg.get_nodeview(c[0]).next() if x.node().node_id()==c[1]]:
            print("Direct connection")
            print(ge.queue_path_detachment([c[0],c[1]],True))
            continue
        pf=SDG.PathFinder(ws,c[0],c[1],9)
        pf.load_lrseqs(conn[1],lrr)
        pf.index_seqs()
        print("seqs loaded")
        pscores=[]
        for p in ws.sdg.find_all_paths_between(c[0],c[1],int(c[2])+1000,1000,False):
            p.nodes=[c[0]]+p.nodes+[c[1]]
            pfsp=SDG.PFScoredPath(pf,c[0],c[1])
            pfsp.path.nodes=p.nodes
            pfsp.find_hits()
            #print("\n",[x for x in p.nodes],pfsp.score(10000))
            #for rh in pfsp.read_hitpos:
                #print("%s"%str(list(rh)))
            #    figure()
            #    plot([x-rh[0] for x in list(rh)])
            #    show()


            s=pfsp.score(100000)
            if s[0]+s[1]: sp=int(10000*s[0]/(s[0]+s[1]))*1.0/100
            else: sp=0
            pscores.append([(sp,s[0],s[1]),p])
        pscores.sort(key=lambda x:x[0],reverse=True)
        for x in pscores[:10]:
            print(x[0],list(x[1].nodes))
            pfsp=SDG.PFScoredPath(pf,c[0],c[1])
            pfsp.path.nodes=x[1].nodes
            pfsp.find_hits()
            #for rh in pfsp.read_hitpos:
            #    #print("%s"%str(list(rh)))
            #    figure()
            #    plot([x for x in list(rh)])
            #    show()
        if pscores:
            #print([c[0]]+pscores[0][1].nodes+[c[1]])
            print(pscores[0][1].nodes)
            #print(ge.queue_path_detachment([c[0]]+pscores[0][1].nodes+[c[1]],True))
            print(ge.queue_path_detachment(pscores[0][1].nodes,True))
        #pf.load_lrseqs(conn[1],lrr)
        #print(pf.lrseqs_as_fasta())
    ge.apply_all()
