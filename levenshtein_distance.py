import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from Levenshtein import distance as lenvenstein_distance


def main():
    train_df = pd.read_csv('novozymes-data/updated_train.csv')

    dist_matrix = np.zeros((len(train_df),len(train_df)))

    for i in tqdm(range(len(train_df))):
        for j in range(max(i-1,0),len(train_df)):

            if abs(len(train_df["protein_sequence"].iloc[i]) - len(train_df["protein_sequence"].iloc[j])) > 15:
                dist_matrix[i,j] = np.inf
            else:
                dist_matrix[i,j] = lenvenstein_distance(train_df["protein_sequence"].iloc[i],
                                                        train_df["protein_sequence"].iloc[j])

            dist_matrix[j,i] = dist_matrix[i,j]

        if i%500 == 0:
            np.save('levenshtein_distance.npy', dist_matrix)

if __name__ == "__main__":
    main()