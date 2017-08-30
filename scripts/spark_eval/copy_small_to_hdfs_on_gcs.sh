#!/usr/bin/env bash

# Copy data prep script to GCS cluster, then run it

gcloud compute scp copy_small_to_hdfs.sh "$GCS_CLUSTER"-m:copy_small_to_hdfs.sh \
  --zone us-central1-a
gcloud compute ssh "$GCS_CLUSTER"-m \
  --command "./copy_small_to_hdfs.sh /user/$USER/small_spark_eval" \
  --zone us-central1-a
