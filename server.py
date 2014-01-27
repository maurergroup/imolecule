"""
A server to communicate with the javascript-based viewer.

Currently implements zeroMQ sockets server side, which are mapped to javascript
websockets wrapped with SocketIO using tornado.
    * zeroMQ - A wrapper around sockets that handles a lot of messiness
               involved with network connections
    * socketIO - A wrapper around javascript websockets that handles the
                 differences in implementations across browser/OS combinations
    * tornado - A Python-based web framework that allows us to convert easily
                between the zeroMQ and socketIO wrappers.
It sounds complicated to use all of these libraries, but it makes this approach
more robust and surprisingly easier.
"""
import os
import multiprocessing
import logging
import traceback

import zmq
from zmq.eventloop import ioloop, zmqstream
ioloop.install()

import tornado
import tornado.web
import tornadio
import tornadio.router
import tornadio.server


class IndexHandler(tornado.web.RequestHandler):

    def get(self):
        self.write(INDEX)


class ClientConnection(tornadio.SocketConnection):

    def on_message(self, message):
        """Evaluates the function pointed to by json-rpc."""
        logging.log(logging.DEBUG, message)
        error = None

        def f(m, q):
            """Wrapper function for worker process to execute."""
            import format_converter
            q.put(getattr(format_converter, m["method"])(**m["params"]))

        # Spawn a thread to protect the server against error signals
        # For larger projects, this could easily be made asynchronous
        queue = multiprocessing.Queue()
        p = multiprocessing.Process(target=f, args=(message, queue))
        p.start()
        p.join(timeout=TIMEOUT)
        if p.exitcode >= 0:
            result = queue.get()
        else:
            error = 1
            if p.exitcode == -1:
                result = traceback.format_exc()
            elif p.exitcode is None:
                result = ("File format conversion timed out! For larger "
                          "files, please run imolecule locally with a large "
                          "`--timeout` parameter.")
            else:
                result = ("Segfault occurred in Open Babel! Please email the "
                          "Open Babel development team with the file used:\n"
                          "<OpenBabel-discuss@lists.sourceforge.net>")
        p.terminate()
        logging.log(logging.DEBUG, result)
        self.send({"result": result, "error": error, "id": message["id"]})


if __name__ == "__main__":
    import argparse
    import webbrowser

    ROOT = os.path.normpath(os.path.dirname(__file__))

    parser = argparse.ArgumentParser(description="Opens a browser-based "
                                     "client that interfaces with the chemical"
                                     " format converter.")
    parser.add_argument("--debug", action="store_true", help="Prints all "
                        "transmitted data streams")
    parser.add_argument("--http-port", type=int, default=8000, help="The port "
                        "on which to serve the website")
    parser.add_argument("--tcp-port", type=int, default=5000, help="The "
                        "server-side tcp connection for python-js interaction")
    parser.add_argument("--timeout", type=int, default=5, help="The maximum "
                        "time, in seconds, allowed for a process to run "
                        "before returning an error")
    args = parser.parse_args()

    HTTP_PORT, TCP_PORT, TIMEOUT = args.http_port, args.tcp_port, args.timeout
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)

    with open("index.html") as index_file:
        INDEX = index_file.read() % {"port": HTTP_PORT}

    WebClientRouter = tornadio.get_router(ClientConnection)
    handler = [(r"/", IndexHandler), WebClientRouter.route(),
               (r'/static/(.*)', tornado.web.StaticFileHandler,
                {'path': ROOT})]
    kwargs = {"enabled_protocols": ["websocket"],
              "socket_io_port": HTTP_PORT}
    application = tornado.web.Application(handler, **kwargs)

    webbrowser.open("http://localhost:%d/" % HTTP_PORT, new=2)

    context = zmq.Context()
    socket = context.socket(zmq.REP)
    socket.bind("tcp://127.0.0.1:%d" % TCP_PORT)
    stream = zmqstream.ZMQStream(socket, tornado.ioloop.IOLoop.instance())
    stream.on_recv(ClientConnection.on_message)
    tornadio.server.SocketServer(application)
