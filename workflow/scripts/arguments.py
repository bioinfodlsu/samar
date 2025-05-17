from typing import List, TypedDict, Union, Optional, Dict, Any, Tuple
from argparse import ArgumentParser


class ArgumentType(TypedDict):
    name: str
    help: str
    type: Optional[float]
    
    
def addArgument(parser_arguments: List[ArgumentType], parser: ArgumentParser) -> ArgumentParser:
  """Adds descriptions to the CLI Arguments

  Args:
      parser_arguments (List[ArgumentType]): List of ArgumentType Dictionaries
      parser (ArgumentParser): The argument parser
  """
  for argument in parser_arguments:
    kwargs = {
        "help": argument['help'],
        "type": argument['type'] if argument['type'] is not None else str
    }
    
    parser.add_argument(argument['name'], **kwargs)
    
  return parser