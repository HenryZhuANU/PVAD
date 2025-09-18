import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_text
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

# Pandas and scikit-learn packages are required to run this file

# Load the CSV file into a pandas DataFrame
data = pd.read_csv('loan_data_set.csv')

# Split the data into input features (X) and target variable (y)
X = data.drop('Loan_Status', axis=1)
y = data['Loan_Status']

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create a decision tree classifier. Use max_depth to control the tree depth
classifier = DecisionTreeClassifier(max_depth=1)

# Train the classifier on the training data
classifier.fit(X_train, y_train)

# Make predictions on the testing data
y_pred = classifier.predict(X_test)

# Print the decision tree model
tree_text = export_text(classifier, feature_names=data.columns[:11].tolist())
print(tree_text)