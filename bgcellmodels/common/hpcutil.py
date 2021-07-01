"""
Tools for working with distributed computing clusters.

@author     Lucas Koelman
@date       3/12/2018
"""

import subprocess
import re

def qstat_get_node_cores(hostname, qstat_line):
    """
    Get node and process IDs from qstat output.
    """
    match = re.search(r"{host}(\d+)/([\d,-]+)".format(host=hostname), qstat_line)
    node_id, process_ids_str = match.groups() # (node-id, process-ids)
    process_ids_int = []
    for processes in process_ids_str.split(','):
        seq = processes.split('-')
        if len(seq) == 0:
            process_ids_int.append(int(seq[0]))
        else:
            process_ids_int.extend(list(range(seq[0], seq[1]+1)))
    return int(node_id), process_ids_int


def qstat_node_sharing_jobs(username, job_id, hostname):
    """
    Get job IDs for jobs running on the same node submitted by user.
    """
    qstat_cmd = "qstat -n -1 -u {user}".format(user=username)
    cmd_output = subprocess.check_output(qstat_cmd, shell=True, universal_newlines=True)
    user_jobs_stats = cmd_output.strip().split()

    # Find node of given job
    job_stats = next((stat for stat in user_jobs_stats if stat.startswith(str(job_id))))
    job_node, job_cores = qstat_get_node_cores(hostname, job_stats)