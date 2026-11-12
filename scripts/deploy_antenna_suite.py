"""Idempotently provision or inspect the dynsyslab AWS Antenna suite.

Run only through the Antenna cloud-operation policy; never place credentials
in configuration files.
"""
from __future__ import annotations
import argparse, json, subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "config" / "antenna_suite.json"
TEMPLATE = ROOT / "infrastructure" / "antenna-suite.yaml"

def aws(config, *args):
    """Run one AWS CLI command with the deployment profile and region."""
    return subprocess.run(["aws", *args, "--profile", "default", "--region", config["region"], "--output", "json"], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

def network(config):
    """Resolve the configured VPC and subnet identifiers for deployment."""
    vpc = json.loads(aws(config, "ec2", "describe-vpcs", "--filters", "Name=isDefault,Values=true", "--query", "Vpcs[0].VpcId").stdout)
    subnet = json.loads(aws(config, "ec2", "describe-subnets", "--filters", f"Name=vpc-id,Values={vpc}", "--query", "sort_by(Subnets,&AvailabilityZone)[0].SubnetId").stdout)
    if not vpc or not subnet: raise RuntimeError("No default VPC/subnet found; provision a network or update this suite design.")
    return vpc, subnet

def main():
    """Deploy the versioned Antenna integration assets after validation."""
    parser=argparse.ArgumentParser(); parser.add_argument("--apply", action="store_true"); parser.add_argument("--status", action="store_true"); args=parser.parse_args()
    config=json.loads(CONFIG.read_text())
    if args.status:
        result=aws(config,"cloudformation","describe-stacks","--stack-name",config["stackName"],"--query","Stacks[0].{Status:StackStatus,Outputs:Outputs}")
        print(result.stdout); return
    if not args.apply: parser.error("--apply is required to provision")
    prior = subprocess.run(["aws","cloudformation","describe-stacks","--stack-name",config["stackName"],"--profile","default","--region",config["region"],"--query","Stacks[0].StackStatus","--output","text"], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if prior.returncode == 0 and prior.stdout.strip() == "ROLLBACK_COMPLETE":
        subprocess.run(["aws","cloudformation","delete-stack","--stack-name",config["stackName"],"--profile","default","--region",config["region"]], check=True)
        subprocess.run(["aws","cloudformation","wait","stack-delete-complete","--stack-name",config["stackName"],"--profile","default","--region",config["region"]], check=True)
    vpc, subnet=network(config)
    command=["aws","cloudformation","deploy","--template-file",str(TEMPLATE),"--stack-name",config["stackName"],"--capabilities","CAPABILITY_NAMED_IAM","--parameter-overrides",f"BucketName={config['bucketName']}",f"VpcId={vpc}",f"SubnetId={subnet}",f"RunnerName={config['runnerName']}",f"InstanceType={config['instanceType']}","--tags","Project=dynsyslab","--profile","default","--region",config["region"]]
    subprocess.run(command,check=True)
    output=json.loads(aws(config,"cloudformation","describe-stacks","--stack-name",config["stackName"],"--query","Stacks[0].Outputs").stdout)
    runner=next(item["OutputValue"] for item in output if item["OutputKey"]=="RunnerInstanceId")
    subprocess.run(["aws","ec2","stop-instances","--instance-ids",runner,"--profile","default","--region",config["region"]],check=True)
    print(json.dumps({"stack":config["stackName"],"bucket":config["bucketName"],"runner":runner,"runnerState":"stopping"}))
if __name__=="__main__": main()
