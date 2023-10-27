#include <iostream>
#include <vector>
#include <random>
#include <stack>
#include <algorithm>
#include <cmath>


class Graph {

public:

    int vertex;
    std::vector<std::vector<int>> matrix;
    std::vector<std::vector<int>> list;

    Graph(int vertex) {
        this->vertex = vertex;
        matrix.resize(vertex, std::vector<int>(vertex, false));
        list.resize(vertex);
    }

    void addVertex(int newVertex) {
        int oldVertex = vertex;


        std::vector<std::vector<int>> newMatrix(newVertex, std::vector<int>(newVertex, false));


        for (int i = 0; i < oldVertex; ++i) {
            for (int j = 0; j < oldVertex; ++j) {
                newMatrix[i][j] = matrix[i][j];
            }
        }

        vertex = newVertex;
        matrix = newMatrix;
    }

    void removeVertex(int v) {
        if (v < 0 || v >= vertex) {
            std::cout << "Invalid vertex index." << std::endl;
            return;
        }

        int oldVertex = vertex;


        std::vector<std::vector<int>> newMatrix(oldVertex - 1, std::vector<int>(oldVertex - 1, false));

        int newRow = 0;
        for (int i = 0; i < oldVertex; ++i) {
            if (i == v) {

                continue;
            }

            int newCol = 0;
            for (int j = 0; j < oldVertex; ++j) {
                if (j == v) {

                    continue;
                }


                newMatrix[newRow][newCol] = matrix[i][j];
                ++newCol;
            }
            ++newRow;
        }

        vertex = oldVertex - 1;
        matrix = newMatrix;
    }


    void printMatrix() {
        for (int i = 0; i < vertex; ++i) {
            for (int j = 0; j < vertex; ++j) {
                std::cout << matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    void printList() {
        std::cout << "Adjacency List:\n";
        for (int i = 0; i < vertex; ++i) {
            std::cout << i << " -> ";
            for (int j : list[i]) {
                std::cout << j << " ";
            }
            std::cout << "\n";
        }
    }

    

};

class undirectedGraph : public Graph {
public:
    undirectedGraph(int vertex) : Graph(vertex) {}

    void addEdge(int source, int destination) {
        if (source >= 0 && source < vertex && destination >= 0 && destination < vertex) {
            matrix[source][destination] = true;
            matrix[destination][source] = true;
        }
    }

    void undirectedGenerator(double edgeProbability) {
        if (edgeProbability < 0.0 || edgeProbability > 1.0) {
            std::cout << "Invalid edge probability." << std::endl;
            return;
        }

        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        for (int i = 0; i < vertex; ++i) {
            for (int j = i + 1; j < vertex; ++j) {
                if (distribution(generator) < edgeProbability) {
                    addEdge(i, j);
                }
            }
        }
    }
    void undirectedtoList() {
        for (int i = 0; i < vertex; ++i) {
            for (int j = 0; j < vertex; ++j)
            {
                if (matrix[i][j] == 1) {
                    list[i].push_back(j);

                }
                else {
                    j = j + 1;
                }
            }
        }
    }


};

class directedGraph : public Graph {
public:
    directedGraph(int vertex) : Graph(vertex) {}

    void addEdge(int source, int destination) {
        if (source >= 0 && source < vertex && destination >= 0 && destination < vertex) {
            matrix[source][destination] = true;
        }
    }

    void directedGenerator(double edgeProbability) {
        if (edgeProbability < 0.0 || edgeProbability > 1.0) {
            std::cout << "Invalid edge probability." << std::endl;
            return;
        }

        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_real_distribution<double> distribution(0.0, 1.0);

        for (int i = 0; i < vertex; ++i) {
            for (int j = 0; j < vertex; ++j) {
                if (i != j && distribution(generator) < edgeProbability) {
                    addEdge(i, j);
                }
            }
        }
    }
    void directedtoList() {
        for (int i = 0; i < vertex; ++i) {
            for (int j = 0; j < vertex; ++j)
            {
                if (matrix[i][j] == 1) {
                    list[i].push_back(j);
                    list[j].push_back(i);

                }
                else {
                    j = j + 1;
                }

            }
        }
    }
};

class weightedGraph : public Graph {
public:
    weightedGraph(int vertex) : Graph(vertex) {}

    void addEdge(int source, int destination, int weight) {
        if (source >= 0 && source < vertex && destination >= 0 && destination < vertex) {
            matrix[source][destination] = weight;
        }
    }

    void weightedGenerator(int minWeight, int maxWeight, double probability) {
        if (minWeight > maxWeight || probability < 0.0 || probability > 1.0) {
            std::cout << "Invalid weight range or probability." << std::endl;
            return;
        }

        std::random_device rd;
        std::default_random_engine generator(rd());
        std::uniform_int_distribution<int> weightDistribution(minWeight, maxWeight);
        std::uniform_real_distribution<double> probDistribution(0.0, 1.0);

        for (int i = 0; i < vertex; ++i) {
            for (int j = 0; j < vertex; ++j) {
                if (i != j) {
                    double randValue = probDistribution(generator);
                    if (randValue <= probability) {
                        int weight = weightDistribution(generator);
                        addEdge(i, j, weight);
                    }
                }
            }
        }
    }

    void weightedtoList() {
        for (int i = 0; i < vertex; ++i) {
            for (int j = 0; j < vertex; ++j)
            {
                if (matrix[i][j] > 1) {
                    list[i].push_back(j);

                }
                else {
                    j = j + 1;
                }

            }
        }
    }

};

void dfs(int v, const std::vector<std::vector<int>>& graph, std::vector<bool>& visited, std::stack<int>& stack) {
    visited[v] = true;

    for (int i = 0; i < graph.size(); ++i) {
        if (graph[v][i] && !visited[i]) {
            dfs(i, graph, visited, stack);
        }
    }

    stack.push(v);
}

std::vector<std::vector<int>> transposeGraph(const std::vector<std::vector<int>>& graph) {
    int numVertices = graph.size();
    std::vector<std::vector<int>> transposed(numVertices, std::vector<int>(numVertices, 0));

    for (int i = 0; i < numVertices; ++i) {
        for (int j = 0; j < numVertices; ++j) {
            transposed[i][j] = graph[j][i];
        }
    }

    return transposed;
}

void dfsUtil(int v, const std::vector<std::vector<int>>& graph, std::vector<bool>& visited, std::vector<int>& component) {
    visited[v] = true;
    component.push_back(v);

    for (int i = 0; i < graph.size(); ++i) {
        if (graph[v][i] && !visited[i]) {
            dfsUtil(i, graph, visited, component);
        }
    }
}

std::vector<std::vector<int>> findStronglyConnectedComponents(const std::vector<std::vector<int>>& graph) {
    int numVertices = graph.size();
    std::vector<bool> visited(numVertices, false);
    std::stack<int> stack;

    for (int i = 0; i < numVertices; ++i) {
        if (!visited[i]) {
            dfs(i, graph, visited, stack);
        }
    }

    std::vector<std::vector<int>> transposed = transposeGraph(graph);

    visited.assign(numVertices, false);
    std::vector<std::vector<int>> stronglyConnectedComponents;

    while (!stack.empty()) {
        int v = stack.top();
        stack.pop();

        if (!visited[v]) {
            std::vector<int> component;
            dfsUtil(v, transposed, visited, component);
            stronglyConnectedComponents.push_back(component);
        }
    }

    return stronglyConnectedComponents;
}



int main() {
    int numVertices = 10;
    double probability = 0.1;
    int numGraphs = 100;

    double totalNumComponents = 0.0;
    double totalAverageSize = 0.0;

    for (int graphIndex = 1; graphIndex <= numGraphs; ++graphIndex) {
        directedGraph graph(numVertices);
        graph.directedGenerator(probability);

        std::vector<std::vector<int>> stronglyConnectedComponents = findStronglyConnectedComponents(graph.matrix);

        int numComponents = stronglyConnectedComponents.size();
        double averageSize = 0.0;
        for (const auto& component : stronglyConnectedComponents) {
            averageSize += component.size();
        }
        if (numComponents > 0) {
            averageSize /= numComponents;
        }

        totalNumComponents += numComponents;
        totalAverageSize += averageSize;

        //  std::cout << "Number of Strongly Connected Components for Graph " << graphIndex << ": " << numComponents << std::endl;
        //  std::cout << "Average Size of Strongly Connected Components for Graph " << graphIndex << ": " << averageSize << std::endl;
    }

    double averageNumComponents = totalNumComponents / numGraphs;
    double averageAverageSize = totalAverageSize / numGraphs;

    std::cout << "Average Number of Strongly Connected Components for " << numGraphs << " Graphs: " << averageNumComponents << std::endl;
    std::cout << "Average Average Size of Strongly Connected Components for " << numGraphs << " Graphs: " << averageAverageSize << std::endl;

    return 0;
}





