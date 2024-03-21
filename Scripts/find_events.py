import argparse
import csv
from collections import defaultdict

def find_repeated_events(file_path):
    event_counts = defaultdict(set)
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            method = row.get('Method')
            start = row.get('Start')
            end = row.get('End')
            if method and start and end:
                event_counts[(start, end)].add(method)

    repeated_events = [(event, methods) for event, methods in event_counts.items() if len(methods) >= 3]
    return repeated_events

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Find repeated events in a CSV file.')
    parser.add_argument('-file', type=str, help='Path to the CSV file containing events data')

    args = parser.parse_args()

    if args.file:
        file_path = args.file
        repeated_events = find_repeated_events(file_path)
        print(repeated_events)
    else:
        print("Please provide the path to the CSV file using the -file argument.")
