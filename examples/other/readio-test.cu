#include "./dataio.h"

int main()
{
    int numberOfNodes = getNumberOfLines("data.txt");
    float4 *nodes = (float4*)malloc((numberOfNodes)*sizeof(float4));;
    loadNodesFromFile("data.txt", nodes);
    
    for(int i = 0; i < numberOfNodes; i++)
    {
        printf("Node %d: (%f, %f)\n", i, nodes[i].x, nodes[i].y);
    }
    
    return(0);
}