def CellPopView(data):
    
    ppr.pprint("CPV: creating data ...")
    adata = createData(data)
    ppr.pprint("CPV: data successfully created ...")
    #ppr.pprint(adata)
    #adata.write("rawData.h5ad")

    #subset by cluster
    Cluster_Key = data['ClusterKey']
    ppr.pprint(Cluster_Key)
    Cluster = data['Cluster']
    Cluster = Cluster
    ppr.pprint(Cluster)
    
    adata = adata[adata.obs[Cluster_Key].isin([Cluster])]
    adata.write("clusterSubset2.h5ad")
    ppr.pprint("CPV: data subset by cluster ...")

    #Split by Condition
    Condition_Key = data['ConditionKey']
    ppr.pprint(Condition_Key)
    Condition1 = data['Cond1']
    ppr.pprint(Condition1)
    Condition2 = data['Cond2']
    ppr.pprint(Condition2)
    adata_1 = adata[adata.obs[Condition_Key].isin([Condition1])]
    adata_2 = adata[adata.obs[Condition_Key].isin([Condition2])]
    ppr.pprint("CPV: data subset by condition ...")
    
    #Extract the data
    Table_1 = adata_1.to_df()
    ppr.pprint(Table_1)
    Table_2 = adata_2.to_df()
    ppr.pprint(Table_2)
    
    Table_1 = Table_1.transpose() #now genes = rows, cells = columns
    Expression_1 = Table_1.mean(axis=1) #extracts the average of every row, e.g. average expression of every gene
    

    Table_2 = Table_2.transpose()
    Expression_2 = Table_2.mean(axis=1) #average expression of every gene in Cluster 3, Condition S

    ppr.pprint("CPV: Expression Values calculated ...")

    #plot and return graph
    #CellPopPlot = 
    #CellPopPlot = plt.figure()
    
    plt.scatter(Expression_1,Expression_2, label = "stars", color = "black", 
            marker = ".",  s =30) 

    #plt.plot([1, 2, 3, 4])

    CellPopPlot = plt.gcf()
    
    #CellPopPlot.savefig("lineplot.png")
    
    ppr.pprint("CPV: Plot created ...")

    #plt.title('Cell Population ' + Cluster) 

    ppr.pprint("CPV: Plot titled...")
    
    #plt.grid()

    ppr.pprint("CPV: Plot gridded...")
    
    #plt.xlabel(Condition1)
    #plt.ylabel(Condition2, rotation = 0)
    
    ppr.pprint("CPV: Plot axis labelled...")

    #CellPopPlot = plt.gcf()
    #CellPopPlot.savefig("CellPop4.png")
    #ppr.pprint(type(CellPopPlot))
    ppr.pprint("CPV: CellPopPlot created...")
    return iostreamFig(CellPopPlot)
