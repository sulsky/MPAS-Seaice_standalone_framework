#-------------------------------------------------------------------------------

def add_output_fields_to_stream(filenameIn,
                                filenameOut,
                                changes,
                                fieldsToAdd=[]):

    try:
        import xml.etree.ElementTree as ET
    except ImportError:
        raise Exception("Module xml.etree.ElementTree needed and not available")

    tree = ET.parse(filenameIn)
    streams = tree.getroot()

    for stream in streams:

        for change in changes:
            if (stream.get('name') == change["streamName"]):
                stream.set(change["attributeName"],change["newValue"])

        for streamsToModify in fieldsToAdd:
            streamName = streamsToModify["streamName"]
            fieldNames = streamsToModify["fieldNames"]
            if (stream.get('name') == streamName):

                for field in list(stream):
                    stream.remove(field)

                for fieldName in fieldNames:
                    field = ET.SubElement(stream, "var")
                    field.set("name",fieldName)

    tree.write(filenameOut)

#-------------------------------------------------------------------------------
