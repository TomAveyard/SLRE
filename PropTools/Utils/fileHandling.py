from io import TextIOWrapper
import os

def openFromRelativePath(relativeFilePath: str) -> TextIOWrapper:

    projectPath = os.path.dirname(os.path.realpath('__file__'))
    absFilePath = os.path.join(projectPath, relativeFilePath)
    
    f = open(absFilePath)

    return f
