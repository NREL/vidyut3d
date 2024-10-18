#!/usr/bin/env bash

MECH_HOME="$(pwd)"
MECH_FILE="${MECH_HOME}/Chemistry.yaml"
bash ../../../../PelePhysics/Support/Mechanism/Models/converter.sh -f "${MECH_FILE}"
