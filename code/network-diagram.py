import matplotlib.pyplot as plt
import networkx as nx

# Define the 12 domains
domains = [
    "Security and Compliance", "Technical Architecture", "Product Delivery",
    "Quality Assurance", "Data Management and Analytics", "Operational Excellence",
    "Team and Culture", "Project and Portfolio Management", "Customer Focus",
    "Innovation and Research", "Vendor and Resource Management", "Performance Measurement"
]

# Create a directed graph
G = nx.DiGraph()

# Add nodes to the graph
for domain in domains:
    G.add_node(domain)

# Define relationships (edges) with labels
relationships = [
    ("Security and Compliance", "Technical Architecture", "informs"),
    ("Technical Architecture", "Product Delivery", "influences"),
    ("Product Delivery", "Quality Assurance", "depends on"),
    # ... other relationships ...
]

# Add edges to the graph
for source, target, label in relationships:
    G.add_edge(source, target, label=label)

# Position nodes using a circular layout
pos = nx.circular_layout(G)

# Draw the graph components
nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='black', node_size=2500, font_size=9)

# Draw edge labels
edge_labels = nx.get_edge_attributes(G, 'label')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=8)

# Show the plot
plt.title('Engineering Maturity Assessment Domains Diagram')
plt.show()
