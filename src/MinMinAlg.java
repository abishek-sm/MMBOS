package ClientPKG;

import ToolPkg.DataLoader;
import ToolPkg.Feat_Normalizer;
import ToolPkg.QSortAlgorithm;
import ToolPkg.SortItem;
import ToolPkg.Statistics;
import ToolPkg.Util;
import java.util.ArrayList;
import java.util.List;


public class MinMinAlg {
    static double corrAVG = 0;
    public static void main(String[] args){
        String filePath, outputFile,scorefile;
        filePath = "C:/Users/cumc/Desktop/dream4size100multifactorial/net5/data.txt";
        String dcfile=       "C:/Users/cumc/Desktop/dream4size100multifactorial/net5/dcormatC.txt";
        outputFile = filePath+".eList.vote.bsdc2475";
        
       double[][] corrMX = DataLoader.loadDoubleData(dcfile); 
        ArrayList<BasicEdge> eList = getRankedEdges_byVoting(corrMX,500000);
    }
    
    public static float[] getEdgesVotes(ArrayList<BasicEdge> candidEdges){
       System.out.println("start getEdgesVotes");
       float[] votes = new float[candidEdges.size()];
       for(int i=0;i<candidEdges.size();i++){
           if(i%1000 ==0)
               System.out.println(i);
           BasicEdge ei = candidEdges.get(i);
           
           for(int j=i+1;j<candidEdges.size();j++){
               BasicEdge ej = candidEdges.get(j);
               int minTriCount = Math.min(ei.allTranglesScores_sortedASCE.length, ej.allTranglesScores_sortedASCE.length);
               double iScore = ei.getPsudMinMinScore(minTriCount);
               double jScore = ej.getPsudMinMinScore(minTriCount);

               if(iScore > jScore)
                   votes[i] +=1;
               else if(jScore > iScore)
                   votes[j] +=1;
               else{
                   votes[j] += 0.5;
                   votes[i] += 0.5;
               }
           }   
       }
        System.out.println("end getEdgesVotes");
       return votes;
    }

    public static void initPsudoMinMinScores(ArrayList<BasicEdge> candidEdges, boolean[][] adjMX, double[][] corrMX){
        System.out.println("start initPsudoMinMinScores");
        for(BasicEdge e:candidEdges){
            double[] triaScores = getEdgeTriangleScores(e.x, e.y, adjMX, corrMX, candidEdges);
            e.setTrianglesScores(triaScores, corrMX[e.x][e.y]);
        }
        System.out.println("end initPsudoMinMinScores");
    }
    public static ArrayList<BasicEdge>  getRankedEdges_byVoting(double[][] corrMX, int initCandidEdges){
        
        ArrayList<BasicEdge> candidEdges =new ArrayList<BasicEdge>();
        double[][] scoreMX = new double[corrMX.length][corrMX.length];
        boolean[][] adjMX = new boolean[corrMX.length][corrMX.length];
        
        initAlg(corrMX, scoreMX, adjMX, candidEdges, initCandidEdges);
        scoreMX =null;
        ArrayList<BasicEdge> edgeSortedList =new ArrayList<BasicEdge>();
        while(candidEdges.size() > 0){
            initPsudoMinMinScores(candidEdges, adjMX, corrMX); 
            float[] voteCount = getEdgesVotes(candidEdges);
            int[] sortedInds = QSortAlgorithm.getIndexesSortASCE(voteCount);
            for(int i=0;i<Math.max(sortedInds.length/1,1);i++){
                BasicEdge e = candidEdges.get(sortedInds[i]);
                edgeSortedList.add(0,e);
                e.depScore = voteCount [sortedInds[i]] ;
                adjMX[e.x][e.y]=false;
                adjMX[e.y][e.x]=false;
            } 
            for(BasicEdge e:edgeSortedList)
                candidEdges.remove(e);
                    
        }
        
        System.out.println("Vote based Ranking Finished");
        return edgeSortedList;   
    }
    
    public static ArrayList<BasicEdge>  getRankedEdges(double[][] corrMX, int initCandidEdges,String scorefile){
        ArrayList<BasicEdge> candidEdges =new ArrayList<BasicEdge>();
        double[][] scoreMX = new double[corrMX.length][corrMX.length];
        boolean[][] adjMX = new boolean[corrMX.length][corrMX.length];
        initAlg(corrMX, scoreMX, adjMX, candidEdges, initCandidEdges);
        
        
        ArrayList<BasicEdge> edgeSortedList =new ArrayList<BasicEdge>();
        System.out.println("Recursive Elimination Start");
        int counter=0;
        
        while(candidEdges.size()>0){
            BasicEdge e = getMinEdge(adjMX, scoreMX, candidEdges);
            e.depScore = scoreMX[e.x][e.y];
            deleteEdge(e, adjMX, scoreMX, corrMX, candidEdges);
            edgeSortedList.add(0,e);
            if((counter++)%1000==0)
                System.out.println("Iteration: "+counter);
        }
        System.out.println("Recursive Elimination Ends");
        DataLoader.file_put_contents(scorefile,scoreMX);
        return edgeSortedList;   
    }
   
    public static void initAlg(double[][] corrMX, double[][] scoreMX, boolean[][] adjMX, ArrayList<BasicEdge> candidEdges, int initCandidEdges){

        corrAVG = Statistics.average(corrMX);
        ArrayList<SortItem> tempCandidList = new ArrayList<SortItem>();
        for(int i=0;i<corrMX.length;i++){
            for(int j=i+1;j<corrMX.length;j++){
                scoreMX[i][j] = corrMX[i][j];
                scoreMX[j][i] = scoreMX[i][j];
                BasicEdge e = new BasicEdge(i, j, scoreMX[i][j]);
                tempCandidList.add(e);
                if(tempCandidList.size()>initCandidEdges*3){
                    QSortAlgorithm.sortDESCE(tempCandidList);
                    while(tempCandidList.size()>initCandidEdges)
                        tempCandidList.remove(tempCandidList.size()-1);
                }
            }   
        }
        
        QSortAlgorithm.sortDESCE(tempCandidList);
        
        for(int i=0;i<initCandidEdges;i++)
            candidEdges.add((BasicEdge)tempCandidList.get(i));
        
        for(int i=0;i<corrMX.length;i++){
            scoreMX[i][i]=-1;
            adjMX[i][i]=false;
            for(int j=i+1;j<corrMX.length;j++){
                scoreMX[i][j]=-1;
                scoreMX[j][i]=-1;
                adjMX[i][j]=false;
                adjMX[j][i]=false;
            }
        }
        for(SortItem si:candidEdges)
        {
            BasicEdge be = (BasicEdge) si;
            adjMX[be.x][be.y] = true;
            adjMX[be.y][be.x] = true;
        }
        
        preComputeMinMinScore(corrMX, scoreMX, adjMX, candidEdges);
    }
 
    public static void preComputeMinMinScore(double[][] corrMX, double[][] scoreMX, boolean[][] adjMX, ArrayList<BasicEdge> candidEdges){
        for(BasicEdge e:candidEdges){
            double score =  getEdgeMinMinScore(e.x, e.y, adjMX, corrMX, candidEdges);
            scoreMX[e.x][e.y] = score;
            scoreMX[e.y][e.x] = score;
        }
    }
    
    public static double[] getEdgeTriangleScores(int x, int y, boolean[][] adjMX, double[][] corrMX, ArrayList<BasicEdge> candidEdges){
        int[] NGs = getCommonNeighbors( x,  y, adjMX);
        if(NGs ==null || NGs.length<1)
            return null;
        
         double[] scores = new double[NGs.length];
         for(int i=0;i<NGs.length;i++){
             scores[i] = corrMX[x][y] - Math.min(corrMX[NGs[i]][x],corrMX[NGs[i]][y]);
         }
         return scores;
    }
    public static double getEdgeMinMinScore(int x, int y, boolean[][] adjMX, double[][] corrMX, ArrayList<BasicEdge> candidEdges){
        double peerwiseScore = corrMX[x][y] ;
        double score = peerwiseScore ;
        double[] trianglesScores = getEdgeTriangleScores( x, y, adjMX, corrMX, candidEdges); 
        if(trianglesScores!=null){
            for(int i=0;i<trianglesScores.length;i++){
                 score = Math.min(score, trianglesScores[i]);
             }    
        }
        return score;
    }
    
    public static int[] getCommonNeighbors(int x, int y, boolean[][] adjMX){
        ArrayList<Integer> list= new ArrayList<Integer>();
        for(int i=0;i<adjMX.length;i++){
            if(i==x || i==y)
                continue;
            
            if(adjMX[i][x]==true && adjMX[i][y]==true)
                list.add(i);
        }
        return Util.toIntArray(list);
    }

    public static void deleteEdge(BasicEdge e, boolean[][] adjMX, double[][] scoreMX, double[][] corrMX, ArrayList<BasicEdge> candidEdges){
        
        adjMX[e.x][e.y]=false;
        adjMX[e.y][e.x]=false;
        BasicEdge.removeEdge(e.x, e.y, candidEdges);
        
        int[] commonNGs = getCommonNeighbors(e.x, e.y, adjMX);
        for(int ind:commonNGs){
            double sc_x = getEdgeMinMinScore(e.x, ind, adjMX, corrMX, candidEdges);
            double sc_y = getEdgeMinMinScore(e.y, ind, adjMX, corrMX, candidEdges);
            scoreMX[ind][e.x] = sc_x;
            scoreMX[e.x][ind] = sc_x;
            scoreMX[ind][e.y] = sc_y;
            scoreMX[e.y][ind] = sc_y;
        }
    }    
    
    static BasicEdge getMinEdge(boolean[][] adjMX, double[][] scoreMX, ArrayList<BasicEdge> candidEdges){
        double minScore = 9999999;
        BasicEdge minE = null;
        for(BasicEdge e: candidEdges){
            if(adjMX[e.x][e.y]==false)
                continue;
            
            if(scoreMX[e.x][e.y]<minScore){
                minScore = scoreMX[e.x][e.y];
                minE =e;
            }
        }
        return minE;
    }
    
    public static void normalizeCorrMX_byMedian(double[][] corr){
        System.out.println("Start CorrMX byMedian Normalization");
        
        double[] medians =new double[corr.length];
        for(int i=0;i<corr.length;i++){
            medians[i] = Statistics.getMedian(corr[i]);
        }
        
        double medianAVG = Statistics.average(medians);
        
        for(int i=0;i<corr.length;i++){
            for(int j=i+1;j<corr.length;j++){
                corr[i][j] = corr[i][j] *  2*medianAVG / (medians[i] + medians[j]);
                corr[j][i] = corr[i][j];
            }  
        }
        
        System.out.println("End CorrMX byMedian Normalization");
    }
    
    public static void normalizeCorrMX_byMedian2(double[][] corr){
        System.out.println("Start CorrMX byMedian Normalization 2");
        
        double[] medians =new double[corr.length];
        for(int i=0;i<corr.length;i++){
            medians[i] = Statistics.getMedian(corr[i]);
        }
        
        double medianAVG = Statistics.average(medians);
        double medianSTD = Statistics.ramiStandardDeviation(medians);
        System.out.println("median(AVG : STD) = ("+medianAVG +" : " + medianSTD+")");
        
        for(int i=0;i<corr.length;i++){
            for(int j=i+1;j<corr.length;j++){
                corr[i][j] = corr[i][j] - medians[i] - medians[j] + 2*medianAVG;
                corr[i][j] = Math.max(0, corr[i][j] );
                corr[j][i] = corr[i][j];
            }  
        }
        
        System.out.println("End CorrMX byMedian Normalization 2");
    }

    public static void normalizeCorrMX_byMean(double[][] corr){
        System.out.println("Start CorrMX byMedian Normalization");
        
        double[] means =new double[corr.length];
        for(int i=0;i<corr.length;i++){
            means[i] = Statistics.average(corr[i]);
        }
        
        double meanAVG = Statistics.average(means);
        
        for(int i=0;i<corr.length;i++){
            for(int j=i+1;j<corr.length;j++){
                //corr[i][j] = corr[i][j] *  meanAVG / Math.sqrt(means[i] * means[j]);
                corr[i][j] = corr[i][j] *  2 * meanAVG / (means[i] + means[j]);
                corr[j][i] = corr[i][j];
            } 
        }
        
        System.out.println("End CorrMX byMedian Normalization");
    }
}

class BasicEdge implements SortItem{
    int x=-1, y=-1;
    double depScore=-1;
    
    double[] allTranglesScores_sortedASCE;
    double[] psudoMinMin_BySort;
    double[] psudoMinMin_ByBootStrap;
    
    BasicEdge(int x, int y, double depScore){
       this.x = Math.min(x, y);
       this.y= Math.max(y,x);
       this.depScore = depScore;
    }
    
    public void setTrianglesScores(double[] allTranglesScores, double basicCorrScore){
        if(allTranglesScores!=null){
            this.allTranglesScores_sortedASCE = new double[allTranglesScores.length+1];
            for(int i=0;i<allTranglesScores.length;i++)
                this.allTranglesScores_sortedASCE[i] = allTranglesScores[i];
            
            this.allTranglesScores_sortedASCE[this.allTranglesScores_sortedASCE.length - 1] = basicCorrScore;
            QSortAlgorithm.sortASCE(this.allTranglesScores_sortedASCE);
        }else{
            this.allTranglesScores_sortedASCE = new double[1];
            this.allTranglesScores_sortedASCE[0] = basicCorrScore;
        }

      preCompute_PsudoMinMin_BySort();
      preCompute_PsudoMinMin_BootStrap();
    }
    
    public void preCompute_PsudoMinMin_BootStrap(){

        psudoMinMin_ByBootStrap = null;
        double[]temp = new double[this.allTranglesScores_sortedASCE.length];
        for(int i=0;i<temp.length;i++){
            if(i>50)
                temp[i] = temp[i-1];
            else
                temp[i] = getPsudMinMinScore_ByBootStrap(i+1);
        }
        psudoMinMin_ByBootStrap = temp;
    }
    
    public void preCompute_PsudoMinMin_BySort(){

        psudoMinMin_BySort = null;
        double[]temp = new double[this.allTranglesScores_sortedASCE.length];
        for(int i=0;i<temp.length;i++)
            temp[i] = getPsudMinMinScore_BySort(i+1);
        
        psudoMinMin_BySort = temp;
    }
    public double getPsudMinMinScore(int maxMinMinCount){
      return getPsudMinMinScore_ByBootStrap(maxMinMinCount);
    }
    
    public double getPsudMinMinScore_ByBootStrap(int maxMinMinCount){
        if(this.psudoMinMin_ByBootStrap!=null)
            return this.psudoMinMin_ByBootStrap[maxMinMinCount-1];
        
        if(maxMinMinCount == allTranglesScores_sortedASCE.length)
            return allTranglesScores_sortedASCE[0];
        
        if(maxMinMinCount == 1)
            return allTranglesScores_sortedASCE[allTranglesScores_sortedASCE.length-1];
        
        maxMinMinCount = Math.min(maxMinMinCount + 0, allTranglesScores_sortedASCE.length);
             
        return Statistics.bootstrappedMin_withReplacement_sortedArr(allTranglesScores_sortedASCE, maxMinMinCount);
       
    }
    
    public double getPsudMinMinScore_BySort(int maxMinMinCount){

        if(maxMinMinCount==1)
          return allTranglesScores_sortedASCE[allTranglesScores_sortedASCE.length-1];
                 
        if(psudoMinMin_BySort!=null)
            return psudoMinMin_BySort[maxMinMinCount-1];
        
        if(maxMinMinCount == allTranglesScores_sortedASCE.length)
            return allTranglesScores_sortedASCE[0];
        
        float maxMinMinCountV = maxMinMinCount;
        float alpha = 1f;
        maxMinMinCountV = alpha*maxMinMinCount + (1-alpha)*allTranglesScores_sortedASCE.length;
            
        //float fracSize = (float)(allTranglesScores_sortedASCE.length-1)/(maxMinMinCount-1);
        float fracSize = (float)(allTranglesScores_sortedASCE.length)/(maxMinMinCountV);
        double scoreSum = 0;
        double scoreProd = 0;
        
        for(int i=0;i<(int)fracSize;i++){
            scoreSum += allTranglesScores_sortedASCE[i];
            if(allTranglesScores_sortedASCE[0] > 0 )
                scoreProd += (1.0/fracSize)*Math.log(allTranglesScores_sortedASCE[i]);
        }
        
        if((fracSize - (int)fracSize) > 0){
            double temp = allTranglesScores_sortedASCE[(int)fracSize];
            scoreSum += temp*(fracSize - (int)fracSize);
            if(allTranglesScores_sortedASCE[0] > 0 )
                scoreProd += (fracSize - (int)fracSize)*Math.log(allTranglesScores_sortedASCE[(int)fracSize]);
            
        }
        return scoreSum / fracSize;
    }
    
    public double getValue(){
        return depScore;
    }
    
    boolean equals(int xv, int yv){
        return (x == Math.min(xv, yv) && y ==Math.max(yv, xv));
    }
    
    @Override
    public boolean equals(Object e0){
        BasicEdge e = (BasicEdge)e0;
        return (x == e.x && y==e.y);
    }
    
    public String toString(){
        return ""+x+"\t"+y+"\t"+depScore;
    }
    
    public static String toString(ArrayList<BasicEdge> eList){
        StringBuffer sb=new StringBuffer("");
        for(BasicEdge e:eList)
            sb.append(e.toString()+"\n");
        
        return sb.toString();
    }
    
    public static void removeEdge(int x, int y, ArrayList<BasicEdge> eList){
        BasicEdge tempE = new BasicEdge(x,y,-1);
        int ind = eList.indexOf(tempE);
        eList.remove(ind);
    }
}
