def CellPopView(data):
    
    adata = createData(data)
    
    #subset by cluster
    adata = adata[adata.obs[data['Cluster]].isin(data['ClusterKey'])]
    
    #Split by Condition
    adata_1 = adata[adata.obs[data['Condition_Key']].isin(data['Cond1'])]
    adata_2 = adata[adata.obs[data['Condition_Key']].isin(data['Cond2'])]
    
    #Extract the data
    Table_1 = adata_1.to_df()
    Table_2 = adata_2.to_df()
    
    Table_1 = Table_1.transpose() #now genes = rows, cells = columns
    Expression_1 = Table_1.mean(axis=1) #extracts the average of every row, e.g. average expression of every gene

    Table_2 = Table_2.transpose()
    Expression_2 = Table_2.mean(axis=1) #average expression of every gene in Cluster 3, Condition S
    
    #plot and return graph
    #CellPopPlot = 
    plt.scatter(Expression_1,Expression_2, label = "stars", color = "black", 
                marker = ".",  s =30) 
    
    plt.title('Cluster ' + str(Cluster)) 
    
    plt.grid()
    
    plt.xlabel(Conditions[0])
    plt.ylabel(Conditions[1], rotation = 0)
    
    CellPopPlot = plt.show()
    return iostreamFig(CellPopPlot)
